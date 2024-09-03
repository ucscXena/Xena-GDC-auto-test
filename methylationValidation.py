import os
import logging
import requests
from concurrent.futures import ThreadPoolExecutor
import json
import subprocess
import pandas
import numpy
import hashlib

logger = logging.getLogger(__name__)


def existing_md5sums(logger, projectName, dataType, sampleDict):
    existingFiles = [file for file in os.listdir(f"gdcFiles/{projectName}/{dataType}") if not file.startswith(".")]
    existingMd5sumFileDict = {md5sum(f"gdcFiles/{projectName}/{dataType}/{file}"): file for file in existingFiles}
    gdcMd5sumFileDict = {sampleDict[sample][fileID]["md5sum"]: fileID for sample in sampleDict for fileID in sampleDict[sample]}
    logger.info(f"{len(gdcMd5sumFileDict)} files found from the GDC for {dataType} data for {projectName}")
    logger.info(f"{len(existingMd5sumFileDict)} files found at gdcFiles/{projectName}/{dataType}")
    fileIdDict = {innerKey: value
                  for outerDict in sampleDict.values()
                  for innerKey, value in outerDict.items()}
    filesNeededToUpdate = {gdcMd5sumFileDict[md5sum]: fileIdDict[gdcMd5sumFileDict[md5sum]]["fileName"] for md5sum in gdcMd5sumFileDict if md5sum not in existingMd5sumFileDict}
    logger.info(f"{len(filesNeededToUpdate)} files needed to update")
    return filesNeededToUpdate
    x = 5


def md5sum(filePath):
    md5_hash = hashlib.md5()

    with open(filePath, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def custom_round(chunk):
    for col in chunk:
        chunk[col] = chunk[col].apply(
            lambda x: numpy.format_float_scientific(x, precision=8) if pandas.notna(x) else numpy.nan)
    
    return chunk


def round_ForNans(x):
    if (pandas.notna(x)):
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
    with open(xenaFile, "r") as xenaData:
        header = xenaData.readline()  # get column labels from xena matrix
        sampleList = header.split("\t")  # split tsv file into list
        sampleList.pop(0)  # remove unnecessary label
        sampleList = [sample.strip() for sample in sampleList]  # make sure there isn't extra whitespace
    
    return sampleList


def getAllSamples(projectName, platform):
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
                        "SeSAMe Methylation Beta Estimation"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        "dna methylation"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [
                        "Methylation Beta Value"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": [
                        "Methylation Array"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.platform",
                    "value": [
                        platform
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


def dataTypeSamples(projectName, samples, platform):
    filesEndpt = "https://api.gdc.cancer.gov/files"
    # MAKE IT SO THAT FILTER GETS DATA TYPE INSERTED
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
                        "SeSAMe Methylation Beta Estimation"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": "dna methylation"
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": [
                        "Methylation Beta Value"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "experimental_strategy",
                    "value": [
                        "Methylation Array"
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
                    "field": "platform",
                    "value": [
                        platform
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
    uniqueSamples = []
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            sampleName = sample["submitter_id"]
            if sampleName in dataTypeDict:
                dataTypeDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                 "md5sum": caseDict["md5sum"]
                                                                 }
            else:
                dataTypeDict[sampleName] = {}
                dataTypeDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                 "md5sum": caseDict["md5sum"]
                                                                 }
                uniqueSamples.append(sampleName)
    
    return dataTypeDict, uniqueSamples


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    splitDF = numpy.array_split(xenaDF.columns, 32)
    tasks = [xenaDF[i] for i in splitDF]
    with ThreadPoolExecutor() as executor:
        result = executor.map(custom_round, tasks)
    # print(executor._max_workers)
    resultDF = pandas.concat(result, axis=1)
    
    return resultDF


def compare(sampleDict, xenaDF, projectName, dataType):
    samplesCorrect = 0
    failed = []
    sampleNum = 1
    total = len(sampleDict)
    for sample in sampleDict:
        fileCount = 0
        for fileID in sampleDict[sample]:
            fileName = sampleDict[sample][fileID]["fileName"]
            sampleFile = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            if fileCount == 0:
                sampleDataDF = pandas.read_csv(sampleFile, sep="\t", names=["compElement", "betaValue"], skiprows=0)
                sampleDataDF.reset_index(inplace=True, drop=True)
                sampleDataDF["nonNanCount"] = 0
                sampleDataDF["nonNanCount"] = sampleDataDF.apply(
                    lambda x: 1 if (not (pandas.isna(x["betaValue"]))) else 0, axis=1)
                sampleDataDF["betaValue"] = sampleDataDF["betaValue"].fillna(0)
            else:
                tempDF = pandas.read_csv(sampleFile, sep="\t", names=["compElement", "betaValue"], skiprows=0)
                tempDF.reset_index(inplace=True, drop=True)
                tempDF["nonNanCount"] = 0
                tempDF["nonNanCount"] = tempDF.apply(lambda x: 1 if (not (pandas.isna(x["betaValue"]))) else 0,
                                                     axis=1)
                tempDF["betaValue"] = tempDF["betaValue"].fillna(0)
                sampleDataDF["nonNanCount"] += tempDF["nonNanCount"]
                sampleDataDF["betaValue"] += tempDF["betaValue"]
            fileCount += 1
        sampleDataDF["betaValue"] = sampleDataDF["betaValue"].astype(float)
        sampleDataDF["betaValue"] = sampleDataDF.apply(lambda x: x["betaValue"]/x["nonNanCount"] if x["nonNanCount"] != 0 else numpy.nan, axis=1)
        sampleDataDF["betaValue"] = sampleDataDF["betaValue"].apply(round_ForNans)
        xenaColumn = xenaDF[sample]
        sampleColumn = sampleDataDF["betaValue"]
        xenaColumn.reset_index(inplace=True, drop=True)
        sampleColumn.reset_index(inplace=True, drop=True)
        if (sampleColumn.equals(xenaColumn)):
            status = "[{:d}/{:d}] Sample: {} - Passed"
            logger.info(status.format(sampleNum, total, sample))
            samplesCorrect += 1
        else:
            status = "[{:d}/{:d}] Sample: {} - Failed"
            logger.info(status.format(sampleNum, total, sample))
            failed.append('{} ({})'.format(sample, sampleNum))
        sampleNum += 1
    
    return failed


def main(projectName, xenaFilePath, dataType):
    logger.info("Testing [{}] data for [{}].".format(dataType, projectName))
    platformDict = {
        "methylation_epic": "illumina methylation epic",
        "methylation27": "illumina human methylation 27",
        "methylation450": "illumina human methylation 450",
        "methylation_epic_v2": "illumina methylation epic v2"
    }
    platform = platformDict[dataType]
    xenaSamples = getXenaSamples(xenaFilePath)
    allSamples = getAllSamples(projectName, platform)
    sampleDict, uniqueSamples = dataTypeSamples(projectName, allSamples, platform)
    xenaDF = xenaDataframe(xenaFilePath)
    if sorted(uniqueSamples) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(uniqueSamples)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in uniqueSamples if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in uniqueSamples]}")
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
    result = compare(sampleDict, xenaDF, projectName, dataType)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return 'PASSED'
    else:
        logger.info("[{}] test failed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return 'FAILED'
