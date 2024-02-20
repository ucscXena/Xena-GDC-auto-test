import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 cnvGeneLevelValidation.py [Project Name] [Xena File Path]")
    exit(1)
projectName = sys.argv[1]
# projectName = "CPTAC-3"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CPTAC-3/Xena_Matrices/CPTAC-3.gene-level_ascat-ngs.tsv"

dataType = "copy_number"
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

    os.system("mkdir gdcFiles")
    os.system("tar -xzf gdcFiles.tar.gz -C gdcFiles")


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
            sampleName = sample["submitter_id"]
            if sample["tissue_type"] == "Tumor":
                if sampleName in dataTypeDict:
                    dataTypeDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]
                else:
                    dataTypeDict[sampleName] = {}
                    dataTypeDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]
    return dataTypeDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    xenaDF = xenaDF.round(10)
    return xenaDF


def compare():
    samplesCorrect = 0
    sampleNum = 0
    for sample in sampleDict:
        fileCount = 0
        print(f"Sample: {sample}\nSample Number: {sampleNum}\n\n")
        for fileID in sampleDict[sample]:
            fileName = sampleDict[sample][fileID]
            sampleFile = "gdcFiles/" + fileID + "/" + fileName
            if fileCount == 0:
                sampleDataDF = pandas.read_csv(sampleFile, sep="\t")
                sampleDataDF["nonNanCount"] = 0
                sampleDataDF["nonNanCount"] = sampleDataDF.apply(lambda x: 1 if (not(pandas.isna(x["copy_number"]))) else 0, axis=1)
                sampleDataDF["copy_number"] = sampleDataDF["copy_number"].fillna(0)
            else:
                tempDF = pandas.read_csv(sampleFile, sep="\t")
                tempDF["nonNanCount"] = 0
                tempDF["nonNanCount"] = tempDF.apply(lambda x: 1 if (not(pandas.isna(x["copy_number"]))) else 0, axis=1)
                tempDF["copy_number"] = tempDF["copy_number"].fillna(0)
                sampleDataDF["nonNanCount"] += tempDF["nonNanCount"]
                sampleDataDF["copy_number"] += tempDF["copy_number"]
            fileCount += 1




        cellsCorrect = 0
        sampleDataDF["copy_number"] = sampleDataDF.apply(lambda x: x["copy_number"]/x["nonNanCount"] if x["nonNanCount"] != 0 else numpy.nan, axis=1)
        sampleDataDF = sampleDataDF.round(10)
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


xenaSamples = getXenaSamples(xenaFilePath)


allSamples = getAllSamples(projectName)
sampleDict = dataTypeSamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)

# print(len(sampleDict), len(xenaSamples))
#
# print(f"Retrieved Samples:\n{sorted(sampleDict)}")
# print(f"Xena Samples:\n{sorted(xenaSamples)}")

if sorted(sampleDict) != sorted(xenaSamples):
    print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples")
    exit(1)

fileIDS = [fileID for sample in sampleDict for fileID in sampleDict[sample]]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")

