import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 geneExpressionValidation.py [Project Name] [Xena File Path]")
    exit(1)

projectName = sys.argv[1]
# projectName = "CPTAC-3"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CPTAC-3/Xena_Matrices/CPTAC-3.star_counts.tsv"

dataType = "unstranded"
workflowType = "STAR - Counts"
dataCategory = "Transcriptome Profiling"
gdcDataType = "Gene Expression Quantification"
experimentalStrategy = "RNA-Seq"


def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    print("Downloading from GDC: ")
    subprocess.run(
        ["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
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
            }
        ]
    }
    params = {
        "filters": json.dumps(dataTypeFilter),
        "fields": "cases.samples.submitter_id,file_id,file_name,cases.samples.tissue_type",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    dataTypeDict = {}
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            sampleName = sample["submitter_id"]
            dataTypeDict[sampleName] = dict(file_id=caseDict["file_id"], file_name=caseDict["file_name"])
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
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=1)
        sampleDataDF = sampleDataDF.drop([0, 1, 2, 3])
        sampleDataDF[dataType] = (numpy.log2(sampleDataDF[dataType] + 1)).round(10)
        for row in range(len(sampleDataDF.index)):
            xenaDataCell = xenaDF.iloc[row][sample]
            sampleDataCell = sampleDataDF.iloc[row][dataType]
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

if sorted(sampleDict) != sorted(xenaSamples):
    print("Samples retrieved from GDC and not in Xena Dataframe")
    print([x for x in sampleDict if x not in xenaSamples])
    print("Samples retrieved from Xena Dataframe and not in GDC retrieved samples")
    print([x for x in xenaSamples if x not in sampleDict])
    exit(1)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")
