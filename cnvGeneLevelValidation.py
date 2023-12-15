import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy

dataType = "copy_number"
projectName = "CGCI-HTMCP-LC"

dataCategory = "copy number variation"
gdcDataType = "Gene Level Copy Number"
experimentalStrategy = "WGS"
workflowType = "AscatNGS"

def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    print("Downloading from GDC: ")
    subprocess.run(["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
                     "Content-Type: application/json", "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"])

    gdcTarfile = tarfile.open("gdcFiles.tar.gz")
    gdcTarfile.extractall("gdcFiles")
    gdcTarfile.close()


def getXenaSamples(xenaFile):  # get all samples from the xena matrix
    with open(xenaFile, "r") as xenaData:
        header = xenaData.readline()  # get column labels from xena matrix
        sampleList = header.split("\t")  # split tsv file into list
        sampleList.pop(0)  # remove unnecessary label
        sampleList = [sample.strip() for sample in sampleList]  # make sure there isn't extra whitespace
    return sampleList


def getAllSamples(projectName):
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
                        dataCategory
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


def dataTypeSamples(samples):
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
                        workflowType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": dataCategory
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
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            if sample["tissue_type"] == "Tumor":
                dataTypeDict[sample["submitter_id"]] = dict(file_id=caseDict["file_id"], file_name=caseDict["file_name"])

    return dataTypeDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    xenaDF = xenaDF.round(10)
    return xenaDF


def compare():
    samplesCorrect = 0
    sampleNum = 0
    for sample in sampleDict:
        print(f"Sample: {sample}\nSample Number: {sampleNum}\n\n")
        cellsCorrect = 0
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t")
        for row in range(len(sampleDataDF.index)):
            xenaDataCell = xenaDF.iloc[row][sample]
            sampleDataCell = sampleDataDF.iloc[row]["copy_number"]
            if (xenaDataCell == sampleDataCell) or (pandas.isna(xenaDataCell) and pandas.isna(sampleDataCell)):
                cellsCorrect += 1
            else:
                print(f"wrong comparison, Sample {sample}, index {row}")
        if cellsCorrect == len(sampleDataDF.index):
            samplesCorrect += 1
        sampleNum += 1
    return samplesCorrect == len(sampleDict)


xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CGCI-HTMCP-LC/Xena_Matrices/CGCI-HTMCP-LC.gene-level_ascat-ngs.tsv"
xenaSamples = getXenaSamples(xenaFilePath)


allSamples = getAllSamples(projectName)
sampleDict = dataTypeSamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)

print(len(sampleDict), len(xenaSamples))

print(f"Retrieved Samples:\n{sorted(sampleDict)}")
print(f"Xena Samples:\n{sorted(xenaSamples)}")
if sorted(sampleDict) != sorted(xenaSamples):
    print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples")
    exit(1)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")