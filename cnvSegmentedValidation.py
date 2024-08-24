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


logger = logging.getLogger(__name__)


def round_ForNans(x):
    if( pandas.notna(x) ):
        return numpy.format_float_scientific(x, precision=10)
    else:
        return numpy.nan


def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    logger.info("Downloading from GDC: ")
    subprocess.run(["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
         "Content-Type: application/json", "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"])

    os.system("mkdir -p gdcFiles")
    os.system("tar -xzf gdcFiles.tar.gz -C gdcFiles")


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
        "fields": "cases.samples.submitter_id,cases.samples.tissue_type,file_id,file_name",
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
                dataTypeDict[sampleName + "." + str(seenDict[sampleName])] = dict(file_id=caseDict["file_id"],
                                                            file_name=caseDict["file_name"])

    return dataTypeDict, list(set(seenSamples))


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t")
    xenaDF["value"] = xenaDF["value"].apply(round_ForNans)
    
    return xenaDF


def sampleDataframe(workflowType, sampleDict):
    dataFrame = pandas.DataFrame()
    # Create a dataframe for all the samples retrieved
    for sample in sampleDict:
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
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
    fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
    downloadFiles(fileIDS)
    # sort data frame
    xenaDF.sort_values(by=sorted(xenaDF), inplace=True)
    # create dataframe for samples
    sampleDf = sampleDataframe(workflowType, sampleDict)
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
