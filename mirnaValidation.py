import requests
import logging
import json
import subprocess
import pandas
import numpy
import os
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
    gdcMd5sumFileDict = {sampleDict[sample][fileID]["md5sum"]: fileID for sample in sampleDict for fileID in
                         sampleDict[sample]}
    logger.info(f"{len(gdcMd5sumFileDict)} files found from the GDC for {dataType} data for {projectName}")
    logger.info(f"{len(existingMd5sumFileDict)} files found at gdcFiles/{projectName}/{dataType}")
    fileIdDict = {innerKey: value
                  for outerDict in sampleDict.values()
                  for innerKey, value in outerDict.items()}
    filesNeededToUpdate = {gdcMd5sumFileDict[md5sum]: fileIdDict[gdcMd5sumFileDict[md5sum]]["fileName"] for md5sum in
                           gdcMd5sumFileDict if md5sum not in existingMd5sumFileDict}
    logger.info(f"{len(filesNeededToUpdate)} files needed to update")
    return filesNeededToUpdate
    x = 5


def round_ForNans(x):
    if( pandas.notna(x) ):
        return numpy.format_float_scientific(x, precision=8)
    else:
        return numpy.nan


'''
Given a list of files a post request will be made to the GDC to download all the files
'''
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


'''
Given a xena matrix file all samples will be extracted in order to compare
if gdc requested samples match
'''
def getXenaSamples(xenaFile):  # get all samples from the xena matrix
    with open(xenaFile, "r") as xenaData:
        header = xenaData.readline()  # get column labels from xena matrix
        sampleList = header.split("\t")  # split tsv file into list
        sampleList.pop(0)  # remove unnecessary label
        sampleList = [sample.strip() for sample in sampleList]  # make sure there isn't extra whitespace
    
    return sampleList


def getAllSamples(projectName, gdcDataType):
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
                        "BCGSC miRNA Profiling"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        "transcriptome profiling"
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
                        "miRNA-Seq"
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(allSamplesFilter),
        "fields": "submitter_sample_ids",
        "format": "json",
        "size": 2000000
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


def miRNASamples(projectName, samples, gdcDataType):
    mirnaSamplesFilter = {
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
                        "BCGSC miRNA Profiling"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": "transcriptome profiling"
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
                        "miRNA-Seq"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.submitter_id",
                    "value": samples
                }
            }
        ]
    }
    filesEndpt = "https://api.gdc.cancer.gov/files"
    params = {
        "filters": json.dumps(mirnaSamplesFilter),
        "fields": "cases.samples.submitter_id,file_id,file_name,md5sum",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    mirnaSamplesDict = {}
    for caseDict in responseJson:
        for submitterDict in caseDict["cases"][0]["samples"]:
            sampleName = submitterDict["submitter_id"]
            if sampleName in mirnaSamplesDict:
                mirnaSamplesDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                     "md5sum": caseDict["md5sum"]
                                                                     }
            else:
                mirnaSamplesDict[sampleName] = {}
                mirnaSamplesDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                     "md5sum": caseDict["md5sum"]
                                                                     }
    
    return mirnaSamplesDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    for column in xenaDF:
        xenaDF[column] = xenaDF[column].apply(round_ForNans)
    return xenaDF

def mirnaDataframe(mirnaSamplesDict, projectName, dataType):
    useColDict = {"mirna": [0, 2],
                  "mirna_isoform": [1, 3]
                  }
    indexColDict = {"mirna": "miRNA_ID",
                  "mirna_isoform": "isoform_coords"
                  }
    useCol = useColDict[dataType]
    indexCol = indexColDict[dataType]
    mirnaDataTitle = "reads_per_million_miRNA_mapped"
    mirnaDataframe = pandas.DataFrame()
    for sample in mirnaSamplesDict:
        sampleDataframe = pandas.DataFrame()
        fileCount = 0
        for fileID in mirnaSamplesDict[sample]:
            fileName = mirnaSamplesDict[sample][fileID]["fileName"]
            sampleFile = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            tempDF = pandas.read_csv(sampleFile, sep="\t", skiprows=0, usecols=useCol, index_col=indexCol)
            tempDF["nonNanCount"] = tempDF.apply(
                lambda x: 1 if (not (pandas.isna(x[mirnaDataTitle]))) else 0, axis=1)
            tempDF[mirnaDataTitle] = tempDF[mirnaDataTitle].fillna(0)
            if fileCount == 0:
                sampleDataframe = tempDF
            else:
                sampleDataframe += tempDF
            fileCount += 1
        sampleDataframe[mirnaDataTitle] = sampleDataframe[mirnaDataTitle].astype(float)
        sampleDataframe[mirnaDataTitle] = sampleDataframe.apply(
            lambda x: x[mirnaDataTitle] / x["nonNanCount"] if x["nonNanCount"] != 0 else numpy.nan, axis=1)
        sampleDataframe[mirnaDataTitle] = numpy.log2(sampleDataframe[mirnaDataTitle] + 1)
        sampleDataframe[mirnaDataTitle] = sampleDataframe[mirnaDataTitle].apply(round_ForNans)
        sampleDataframe.drop("nonNanCount", inplace=True, axis=1)
        sampleDataframe.rename(columns={"reads_per_million_miRNA_mapped": sample}, inplace=True)
        mirnaDataframe = pandas.concat([mirnaDataframe, sampleDataframe], axis=1)
    return mirnaDataframe


def compare(logger, gdcDF, xenaDF):
    failed = []
    sampleNum = 1
    total = len(gdcDF.columns)
    gdcDF.sort_index(inplace=True)
    xenaDF.sort_index(inplace=True)
    for sample in list(gdcDF.columns):
        xenaColumn = xenaDF[sample]
        gdcColumn = gdcDF[sample]
        if not xenaColumn.equals(gdcColumn):
            status = "[{:d}/{:d}] Sample: {} - Failed"
            logger.info(status.format(sampleNum, total, sample))
            failed.append('{} ({})'.format(sample, sampleNum))
        else:
            status = "[{:d}/{:d}] Sample: {} - Passed"
            logger.info(status.format(sampleNum, total, sample))
        sampleNum += 1
    return failed


def main(projectName, xenaFilePath, dataType):
    gdcDataTypeDict = {"mirna": "miRNA Expression Quantification",
                       "mirna_isoform": "Isoform Expression Quantification"}
    gdcDataType = gdcDataTypeDict[dataType]
    logger.info("Testing [{}] data for [{}].".format(dataType, projectName))
    xenaSamples = getXenaSamples(xenaFilePath)
    allSamples = getAllSamples(projectName, gdcDataType)
    mirnaSamplesDict = miRNASamples(projectName, allSamples, gdcDataType)
    xenaDF = xenaDataframe(xenaFilePath)
    if sorted(mirnaSamplesDict) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(mirnaSamplesDict)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in mirnaSamplesDict if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in mirnaSamplesDict]}")
        exit(1)
    if os.path.isdir(f"gdcFiles/{projectName}/{dataType}"):
        fileIDs = existing_md5sums(logger, projectName, dataType, mirnaSamplesDict)
    else:
        fileIDs = [fileID for sample in mirnaSamplesDict for fileID in mirnaSamplesDict[sample]]
        logger.info(f"{len(fileIDs)} files found from the GDC for {dataType} data for {projectName}")
        logger.info(f"0 files found at gdcFiles/{projectName}/{dataType}")
        logger.info(f"{len(fileIDs)} files needed to download")
    if len(fileIDs) != 0:
        downloadFiles(fileIDs, projectName, dataType)
    gdcDF = mirnaDataframe(mirnaSamplesDict, projectName, dataType)
    result = compare(logger, gdcDF, xenaDF)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return 'PASSED'
    else:
        logger.info("[{}] test failed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return 'FAILED'
