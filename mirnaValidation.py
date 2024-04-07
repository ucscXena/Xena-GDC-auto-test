import requests
import json
import tarfile
import subprocess
import pandas
import numpy
import sys
import os
from math import log10, floor


if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 mirnaValidation.py [Project Name] [Xena File Path]")
    exit(1)

projectName = sys.argv[1]
# projectName = "CPTAC-3"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CPTAC-3/Xena_Matrices/CPTAC-3.mirna.tsv"

workflowType = "BCGSC miRNA Profiling"
dataCategory = "Transcriptome Profiling"
gdcDataType = "miRNA Expression Quantification"
experimentalStrategy = "miRNA-Seq"

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

'''
Given a list of files a post request will be made to the GDC to download all the files
'''
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

'''

'''
def getAllSamples():
    casesEndpt = "https://api.gdc.cancer.gov/cases"

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


def miRNASamples(samples):
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
    filesEndpt = "https://api.gdc.cancer.gov/files"
    params = {
        "filters": json.dumps(mirnaSamplesFilter),
        "fields": "cases.samples.submitter_id,file_id,file_name",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    mirnaSamplesDict = {}
    for caseDict in responseJson:
        for submitterDict in caseDict["cases"][0]["samples"]:
            sampleName = submitterDict["submitter_id"]
            mirnaSamplesDict[sampleName] = dict(file_id=caseDict["file_id"], file_name=caseDict["file_name"])
    return mirnaSamplesDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    for column in xenaDF:
        xenaDF[column] = xenaDF[column].apply(round_ForNans)
    return xenaDF

def compare():
    samplesCorrect = 0
    mirnaDataTitle = "reads_per_million_miRNA_mapped"
    for sample in mirnaSamplesDict:
        cellsCorrect = 0
        fileId = mirnaSamplesDict[sample]["file_id"]
        fileName = mirnaSamplesDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", index_col=0)
        sampleDataDF[mirnaDataTitle] = (numpy.log2(sampleDataDF[mirnaDataTitle] + 1)).apply(round_ForNans)
        for row in range(len(sampleDataDF.index)):
            xenaDataCell = xenaDF.iloc[row][sample]
            sampleDataCell = sampleDataDF.iloc[row][mirnaDataTitle]
            if (xenaDataCell == sampleDataCell) or (pandas.isna(xenaDataCell) and pandas.isna(sampleDataCell)):
                cellsCorrect += 1
            else:
                print(f"wrong comparison, Sample {sample}, index {row}")
        if cellsCorrect == len(sampleDataDF.index):
            samplesCorrect += 1
    return samplesCorrect == len(mirnaSamplesDict)


xenaSamples = getXenaSamples(xenaFilePath)

allSamples = getAllSamples()
mirnaSamplesDict = miRNASamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)


if sorted(mirnaSamplesDict) != sorted(xenaSamples):
    print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples")
    print(f"Length of allSamples: {len(allSamples)}")
    print(f"Length of xena samples: {len(xenaSamples)}")
    print(f"Length of mirna samples: {len(mirnaSamplesDict)}")

    print(f"Samples in xena samples and not in mirnaSamplesDict:\n{[x for x in xenaSamples if x not in mirnaSamplesDict]}")
    print(f"Samples in mirnaSamplesDict and not in xenaSamples:\n{[x for x in mirnaSamplesDict if x not in xenaSamples]}")

    exit(1)

fileIDS = [mirnaSamplesDict[x]["file_id"] for x in mirnaSamplesDict]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")


