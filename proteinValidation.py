import pandas
import requests
import logging
import subprocess
import os
import numpy
import hashlib


logger = logging.getLogger(__name__)


def md5sum(filePath):
    md5_hash = hashlib.md5()

    with open(filePath, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def existing_md5sums(logger, projectName, dataType, proteinSamplesDict):
    existingFiles = [file for file in os.listdir(f"gdcFiles/{projectName}/{dataType}") if not file.startswith(".")]
    existingMd5sumFileDict = {md5sum(f"gdcFiles/{projectName}/{dataType}/{file}"): file for file in existingFiles}
    gdcMd5sumFileDict = {proteinSamplesDict[sample][fileID]["md5sum"]: fileID for sample in proteinSamplesDict for fileID in proteinSamplesDict[sample]}
    logger.info(f"{len(gdcMd5sumFileDict)} files found from the GDC for {dataType} data for {projectName}")
    logger.info(f"{len(existingMd5sumFileDict)} files found at gdcFiles/{projectName}/{dataType}")
    fileIdDict = {innerKey: value
                  for outerDict in proteinSamplesDict.values()
                  for innerKey, value in outerDict.items()}
    filesNeededToUpdate = {gdcMd5sumFileDict[md5sum]: fileIdDict[gdcMd5sumFileDict[md5sum]]["fileName"] for md5sum in gdcMd5sumFileDict if md5sum not in existingMd5sumFileDict}
    logger.info(f"{len(filesNeededToUpdate)} files needed to update")
    return filesNeededToUpdate
    x = 5


def getXenaSamples(xenaFile):
    with open(xenaFile, "r") as xenaData:
        header = xenaData.readline()  # get column labels from xena matrix
        sampleList = header.split("\t")  # split tsv file into list
        sampleList.pop(0)  # remove unnecessary label
        sampleList = [sample.strip() for sample in sampleList]  # make sure there isn't extra whitespace

    return sampleList


def proteinSamples(projectName):
    filesEndpoint = "https://api.gdc.cancer.gov/files"
    proteinSamplesFilter = {
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
                    "field": "data_category",
                    "value": [
                        "proteome profiling"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": [
                        "Protein Expression Quantification"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "experimental_strategy",
                    "value": [
                        "Reverse Phase Protein Array"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "platform",
                    "value": [
                        "rppa"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "access",
                    "value": [
                        "open"
                    ]
                }
            }
        ]
    }
    params = {
        "filters": proteinSamplesFilter,
        "fields": "cases.samples.submitter_id,file_id,file_name,md5sum",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpoint, json=params, headers={"Content-Type": "application/json"})
    responseJson = response.json()["data"]["hits"]
    proteinSamplesDict = {}
    for caseDict in responseJson:
        for submitterDict in caseDict["cases"][0]["samples"]:
            sampleName = submitterDict["submitter_id"]
            if sampleName in proteinSamplesDict:
                proteinSamplesDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                       "md5sum": caseDict["md5sum"]
                                                                       }
            else:
                proteinSamplesDict[sampleName] = {}
                proteinSamplesDict[sampleName][caseDict["file_id"]] = {"fileName": caseDict["file_name"],
                                                                       "md5sum": caseDict["md5sum"]
                                                                       }

    return proteinSamplesDict


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


def proteinDataframe(proteinSamples, projectName, dataType):
    proteinDataframe = pandas.DataFrame()
    for sample in proteinSamples:
        sampleDataframe = pandas.DataFrame()
        fileCount = 0
        for fileName in [x["fileName"] for x in list(proteinSamples[sample].values())]:
            filePath = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            tempDF = pandas.read_csv(filePath, sep="\t", usecols=[4, 5], index_col=0)
            tempDF["nonNanCount"] = tempDF.apply(
                    lambda x: 1 if (not (pandas.isna(x["protein_expression"]))) else 0, axis=1)
            tempDF.fillna(0, inplace=True)
            fileCount += 1
            if fileCount == 1:
                sampleDataframe = tempDF
            else:
                sampleDataframe += tempDF
        sampleDataframe["protein_expression"] = sampleDataframe.apply(
            lambda x: x["protein_expression"] / x["nonNanCount"] if x["nonNanCount"] != 0 else numpy.nan, axis=1)
        sampleDataframe.drop("nonNanCount", inplace=True, axis=1)
        sampleDataframe.rename(columns={"protein_expression": sample}, inplace=True)
        proteinDataframe[sample] = sampleDataframe[sample]
    return proteinDataframe


def compare(logger, gdcDF, xenaDF):
    failed = []
    sampleNum = 1
    total = len(gdcDF.columns)
    for sample in gdcDF:
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
    xenaSamples = getXenaSamples(xenaFilePath)
    proteinSamplesDict = proteinSamples(projectName)
    if sorted(proteinSamplesDict) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(proteinSamplesDict)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in proteinSamplesDict if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in proteinSamplesDict]}")
        exit(1)
    if os.path.isdir(f"gdcFiles/{projectName}/{dataType}"):
        fileIDs = existing_md5sums(logger, projectName, dataType, proteinSamplesDict)
    else:
        fileIDs = [fileID for sample in proteinSamplesDict for fileID in proteinSamplesDict[sample]]
        logger.info(f"{len(fileIDs)} files found from the GDC for {dataType} data for {projectName}")
        logger.info(f"0 files found at gdcFiles/{projectName}/{dataType}")
        logger.info(f"{len(fileIDs)} files needed to download")
    if len(fileIDs) != 0:
        downloadFiles(fileIDs, projectName, dataType)
    xenaDF = pandas.read_csv(xenaFilePath, sep="\t", index_col=0)
    gdcDF = proteinDataframe(proteinSamplesDict, projectName, dataType)

    logger.info("Testing in progress ...")

    result = compare(logger, gdcDF, xenaDF)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return "PASSED"
    else:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return "FAILED"

