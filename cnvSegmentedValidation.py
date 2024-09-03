import os
import sys
import logging
import pandas
import difflib
import requests
import json
import subprocess
import pandas
import numpy
import hashlib


logger = logging.getLogger(__name__)


def md5sum(filePath):
    md5_hash = hashlib.md5()

    with open(filePath, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def existing_md5sums(logger, projectName, dataType, sampleDict):
    existingFiles = [file for file in os.listdir(f"gdcFiles/{projectName}/{dataType}") if not file.startswith(".")]
    existingMd5sumFileDict = {md5sum(f"gdcFiles/{projectName}/{dataType}/{file}"): file for file in existingFiles}
    gdcMd5sumFileDict = {sampleDict[sample][fileID]["md5sum"]: fileID for sample in sampleDict for fileID in sampleDict[sample]}
    logger.info(f"{len(gdcMd5sumFileDict)} files found from the GDC for {dataType} data for {projectName}")
    logger.info(f"{len(existingMd5sumFileDict)} files found at gdcFiles/{projectName}/{dataType}")
    fileIdDict = {innerKey: value
                  for outerDict in sampleDict.values()
                  for innerKey, value in outerDict.items()}
    filesNeededToUpdate = {gdcMd5sumFileDict[md5sum]: fileIdDict[gdcMd5sumFileDict[md5sum]]["file_name"] for md5sum in gdcMd5sumFileDict if md5sum not in existingMd5sumFileDict}
    logger.info(f"{len(filesNeededToUpdate)} files needed to update")
    return filesNeededToUpdate


def round_ForNans(x):
    if( pandas.notna(x) ):
        return numpy.format_float_scientific(x, precision=8)
    else:
        return numpy.nan


def downloadFiles(fileList, projectName, dataType):
    if isinstance(fileList, list):
        ids = fileList
    elif isinstance(fileList, dict):
        ids = list(fileList.keys())
    jsonPayload = {
        "ids": ids
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    logger.info("Downloading from GDC: ")

    outputDir = f"gdcFiles/{projectName}/{dataType}"
    os.makedirs(outputDir, exist_ok=True)

    curlCommand = [
        "curl", "--request", "POST", "--header", "Content-Type: application/json",
        "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"
    ]

    if len(fileList) != 1:
        outputFile = "gdcFiles.tar.gz"
        curlCommand.extend(["-o", outputFile])
        subprocess.run(curlCommand)
        os.system(f"tar --strip-components=1 -xzf  gdcFiles.tar.gz -C {outputDir}")
    else:
        outputFile = f"{outputDir}/{list(fileList.values())[0]}"
        curlCommand.extend(["-o", outputFile])
        subprocess.run(curlCommand)


def getXenaSamples(xenaFile):  # get all samples from the xena matrix
    xenaMatrix = pandas.read_csv(xenaFile, sep="\t")
    sampleList = list(xenaMatrix["sample"].unique())
    
    return sampleList


def getAllSamples(projectName, workflowType, gdcDataType, experimentalStrategy):
    casesEndpt = "https://api.gdc.cancer.gov/cases"
    allSamplesFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        projectName
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": [
                        workflowType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        "copy number variation"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [
                        gdcDataType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": [
                        experimentalStrategy
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(allSamplesFilter),
        "fields": "submitter_sample_ids",
        "format": "json",
        "size": 20000
    }
    response = requests.post(casesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    allSamples = []
    for caseDict in responseJson:
        for sample in caseDict["submitter_sample_ids"]:
            allSamples.append(sample)
    
    return allSamples


def unpeelJson(jsonObj):
    jsonObj = jsonObj.get("data").get("hits")
    
    return jsonObj


def dataTypeSamples(projectName, workflowType, gdcDataType, experimentalStrategy, samples):
    filesEndpt = "https://api.gdc.cancer.gov/files"
    dataTypeFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        projectName
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": [
                        workflowType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": "copy number variation"
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": [
                        gdcDataType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "experimental_strategy",
                    "value": [
                        experimentalStrategy
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.submitter_id",
                    "value": samples
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.tissue_type",
                    "value": [
                        "tumor"
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(dataTypeFilter),
        "fields": "cases.samples.submitter_id,cases.samples.tissue_type,file_id,file_name,md5sum",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    dataTypeDict = {}
    # create seen dict to see how many times a sample has been seen
    seenDict = {}
    seenSamples = []
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            sampleName = sample["submitter_id"]
            if sample["tissue_type"] == "Tumor":
                seenDict[sampleName] = seenDict.get(sampleName, 0) + 1
                seenSamples.append(sampleName)
                dataTypeDict[sampleName + "." + str(seenDict[sampleName])] = {caseDict["file_id"] : {"file_name": caseDict["file_name"],
                                                                                                                "md5sum": caseDict["md5sum"]}}

    return dataTypeDict, list(set(seenSamples))


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t")
    xenaDF["value"] = xenaDF["value"].apply(round_ForNans)
    
    return xenaDF


def sampleDataframe(workflowType, sampleDict, projectName, dataType):
    dataFrame = pandas.DataFrame()
    # Create a dataframe for all the samples retrieved
    for sample in sampleDict:
        for fileID in sampleDict[sample]:
            fileName = sampleDict[sample][fileID]["file_name"]
            sampleFile = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            normalSampleName = sample[:sample.index(".")]
            # Create data frame for sample data
            sampleDataDF = pandas.read_csv(sampleFile, sep="\t")
            sampleDataDF.rename(columns={'Chromosome': 'Chrom'}, inplace=True)
            sampleDataDF.rename(columns={'GDC_Aliquot': 'sample'}, inplace=True)
            if( workflowType == "DNAcopy" ):
                sampleDataDF.rename(columns={'Segment_Mean': 'value'}, inplace=True)
            elif( workflowType == "AscatNGS" or workflowType == "ASCAT2" or workflowType == "ASCAT3"):
                sampleDataDF.rename(columns={'Copy_Number': 'value'}, inplace=True)
            sampleDataDF.drop(columns=['Major_Copy_Number', 'Minor_Copy_Number', 'Num_Probes'], inplace=True, errors="ignore")
            sampleDataDF.replace(sampleDataDF.iloc[0].iat[0], normalSampleName, inplace=True)
            dataFrame = pandas.concat([dataFrame, sampleDataDF])
    dataFrame["value"] = dataFrame["value"].apply(round_ForNans)
    
    return dataFrame


def main(projectName, xenaFilePath, dataType):
    logger.info("Testing [{}] data for [{}].".format(dataType, projectName))
    workflowDict = {
        "masked_cnv_DNAcopy": "DNAcopy",
        "segment_cnv_ascat-ngs": "AscatNGS",
        "allele_cnv_ascat2": "ASCAT2",
        "allele_cnv_ascat3": "ASCAT3",
        "segment_cnv_DNAcopy": "DNAcopy"
    }
    gdcDataTypeDict = {
        "masked_cnv_DNAcopy": "Masked Copy Number Segment",
        "segment_cnv_ascat-ngs": "Copy Number Segment",
        "allele_cnv_ascat2": "Allele-specific Copy Number Segment",
        "allele_cnv_ascat3": "Allele-specific Copy Number Segment",
        "segment_cnv_DNAcopy": "Copy Number Segment"
    }
    experimentalStrategyDict = {
        "masked_cnv_DNAcopy": "Genotyping Array",
        "segment_cnv_ascat-ngs": "WGS",
        "allele_cnv_ascat2": "Genotyping Array",
        "allele_cnv_ascat3": "Genotyping Array",
        "segment_cnv_DNAcopy": "Genotyping Array"
    }
    workflowType = workflowDict[dataType]
    gdcDataType = gdcDataTypeDict[dataType]
    experimentalStrategy = experimentalStrategyDict[dataType]
    xenaSamples = getXenaSamples(xenaFilePath)
    allSamples = getAllSamples(projectName, workflowType, gdcDataType, experimentalStrategy)
    sampleDict, seenSamples = dataTypeSamples(projectName, workflowType, gdcDataType, experimentalStrategy, allSamples)
    xenaDF = xenaDataframe(xenaFilePath)
    if sorted(seenSamples) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(seenSamples)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in seenSamples if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in seenSamples]}")
        exit(1)
    if os.path.isdir(f"gdcFiles/{projectName}/{dataType}"):
        fileIDs = existing_md5sums(logger, projectName, dataType, sampleDict)
    else:
        fileIDs = [fileID for sample in sampleDict for fileID in sampleDict[sample]]
        logger.info(f"{len(fileIDs)} files found from the GDC for {dataType} data for {projectName}")
        logger.info(f"0 files found at gdcFiles/{projectName}/{dataType}")
        logger.info(f"{len(fileIDs)} files needed to download")
    if len(fileIDs) != 0:
        downloadFiles(fileIDs, projectName, dataType)
    # sort data frame
    xenaDF.sort_values(by=sorted(xenaDF), inplace=True)
    # create dataframe for samples
    sampleDf = sampleDataframe(workflowType, sampleDict, projectName, dataType)
    # sort sample dataframe as well
    sampleDf.sort_values(by=sorted(sampleDf), inplace=True)
    # then reset index ordering for each one
    xenaDF.reset_index(inplace=True, drop=True)
    sampleDf.reset_index(inplace=True, drop=True)
    with open("sampleDF.csv", "w") as sampleFile:
        sampleDf.to_csv(sampleFile)
    with open("xenaDF.csv", "w") as xenaDfFile:
        xenaDF.to_csv(xenaDfFile)
    if sampleDf.equals(xenaDF):
        logger.info("Testing in progress ...")
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return 'PASSED'
    else:
        logger.info("[{}] test failed for [{}].".format(dataType, projectName))
        logger.info("Diff file is being generated with unequal values.")
        with open("sampleDF.csv", "r") as sampleFile:
            with open("xenaDF.csv", "r") as xenaDfFile:
                # if they are not equal then output diff of both files
                sys.stdout.writelines(difflib.unified_diff(sampleFile.readlines(), xenaDfFile.readlines(),
                                                        fromfile="sampleDF.csv", tofile="xenaDF.csv"))                                        
        return 'FAILED'
