import pandas
import requests
import logging
import subprocess
import os
import numpy


logger = logging.getLogger(__name__)


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
        "fields": "cases.samples.submitter_id,file_id,file_name",
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
                proteinSamplesDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]
            else:
                proteinSamplesDict[sampleName] = {}
                proteinSamplesDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]

    return proteinSamplesDict


def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))
    logger.info("Downloading from GDC: ")
    subprocess.run(["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST",
                    "--header", "Content-Type: application/json", "--data", "@payload.txt",
                    "https://api.gdc.cancer.gov/data"])
    os.system("mkdir -p gdcFiles")
    os.system("tar -xzf gdcFiles.tar.gz -C gdcFiles")


def proteinDataframe(proteinSamples):
    proteinDataframe = pandas.DataFrame()
    for sample in proteinSamples:
        sampleDataframe = pandas.DataFrame()
        fileCount = 0
        for fileID, fileName in proteinSamples[sample].items():
            filePath = "gdcFiles/{}/{}".format(fileID, fileName)
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
    fileIDs = [fileID for sample in proteinSamplesDict for fileID in proteinSamplesDict[sample]]
    # downloadFiles(fileIDs)
    xenaDF = pandas.read_csv(xenaFilePath, sep="\t", index_col=0)
    gdcDF = proteinDataframe(proteinSamplesDict)

    logger.info("Testing in progress ...")

    result = compare(logger, gdcDF, xenaDF)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return "PASSED"
    else:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return "FAILED"

