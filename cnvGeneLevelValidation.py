import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy
from math import log10, floor

if ( len(sys.argv) != 4 ):
    print("Usage:\npython3 cnvGeneLevelValidation.py [Project Name] [Xena File Path] [Workflow Type]")
    print("Valid Workflow Type: ['ABSOLUTE LiftOver', 'ASCAT2', 'ASCAT3', 'ASCATNGS']")
    exit(1)
projectName = sys.argv[1]
# projectName = "CGCI-HTMCP-LC"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CGCI-HTMCP-LC/Xena_Matrices/CGCI-HTMCP-LC.gene-level_ascat-ngs.tsv"
workflowType = sys.argv[3]

experimentalStrategyDict = {
    "ABSOLUTE LiftOver": "Genotyping Array",
    "ASCAT2": "Genotyping Array",
    "ASCAT3": "Genotyping Array",
    "ASCATNGS": "WGS"
}

if workflowType not in experimentalStrategyDict:
    print("Invalid Workflow Type")
    print("Valid Workflow Types: ['ABSOLUTE LiftOver', 'ASCAT2', 'ASCAT3', 'ASCATNGS']")
    exit(1)

dataType = "copy_number"
dataCategory = "Copy Number Variation"
gdcDataType = "Gene Level Copy Number"
experimentalStrategy = experimentalStrategyDict[workFlowType]

def round_ForNans(x):
    if( pandas.notna(x) ):
        return numpy.format_float_scientific(x, precision=10)
    else:
        return numpy.nan

# From https://github.com/corriander/python-sigfig/blob/dev/sigfig/sigfig.py
def round_(x, n):
    """Round a float, x, to n significant figures.

	Caution should be applied when performing this operation.
	Significant figures are an implication of precision; arbitrarily
	truncating floats mid-calculation is probably not Good Practice in
	almost all cases.

	Rounding off a float to n s.f. results in a float. Floats are, in
	general, approximations of decimal numbers. The point here is that
	it is very possible to end up with an inexact number:

		roundsf(0.0012395, 3)
		0.00124
	    roundsf(0.0012315, 3)
		0.0012300000000000002

	Basically, rounding in this way probably doesn't do what you want
	it to.
    """
    n = int(n)
    x = float(x)

    if x == 0: return 0

    e = floor(log10(abs(x)) - n + 1)  # exponent, 10 ** e
    shifted_dp = x / (10 ** e)  # decimal place shifted n d.p.
    return round(shifted_dp) * (10 ** e)  # round and revert

def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    print("Downloading from GDC: ")
    subprocess.run(["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
                     "Content-Type: application/json", "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"])

    os.system("mkdir -p gdcFiles")
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
    for column in xenaDF:
        if pandas.api.types.is_numeric_dtype(xenaDF[column]):
            xenaDF[column] = xenaDF[column].apply(round_ForNans)
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
        sampleDataDF["copy_number"] = sampleDataDF["copy_number"].apply(round_ForNans)
        xenaColumn = xenaDF[sample]
        sampleColumn = sampleDataDF["copy_number"]
        xenaColumn.reset_index(inplace=True, drop=True)
        sampleColumn.reset_index(inplace=True, drop=True)
        if (sampleColumn.equals(xenaColumn)):
            print("success")
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
    print("Samples retrieved from GDC and not in Xena Dataframe")
    print([x for x in sampleDict if x not in xenaSamples])
    print("Samples retrieved from Xena Dataframe and not in GDC retrieved samples")
    print([x for x in xenaSamples if x not in sampleDict])
    exit(1)

fileIDS = [fileID for sample in sampleDict for fileID in sampleDict[sample]]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")

