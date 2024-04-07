import os
import sys

import pandas as pd
import difflib
import requests
import json
import subprocess
import tarfile
import pandas
import numpy
from math import floor, log10

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 cnvSegmentedValidation.py [Project Name] [Xena File Path]")
    exit(1)
projectName = sys.argv[1]
# projectName = "HCMI-CMDC"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/HCMI-CMDC/Xena_Matrices/HCMI-CMDC.segment_cnv_ascat-ngs.tsv"

dataCategory = "copy number variation"
gdcDataType = "Copy Number Segment"
experimentalStrategy = "WGS"
workflowType = "AscatNGS"

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
    xenaMatrix = pandas.read_csv(xenaFile, sep="\t")
    sampleList = list(xenaMatrix["sample"].unique())
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



def sampleDataframe():
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
        sampleDataDF.rename(columns={'Copy_Number': 'value'}, inplace=True)
        sampleDataDF.drop(columns=['Major_Copy_Number', 'Minor_Copy_Number'], inplace=True)
        sampleDataDF.replace(sampleDataDF.iloc[0].iat[0], normalSampleName, inplace=True)
        dataFrame = pandas.concat([dataFrame, sampleDataDF])
    dataFrame["value"] = dataFrame["value"].apply(round_ForNans)
    return dataFrame


xenaSamples = getXenaSamples(xenaFilePath)

allSamples = getAllSamples(projectName)
sampleDict, seenSamples = dataTypeSamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)

# print(len(sampleDict), len(xenaSamples))
#
# print(f"Retrieved Samples:\n{sorted(sampleDict)}")
# print(f"Xena Samples:\n{sorted(xenaSamples)}")
if sorted(seenSamples) != sorted(xenaSamples):
    print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples")
    exit(1)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)
# sort data frame
xenaDF.sort_values(by=sorted(xenaDF), inplace=True)

# create dataframe for samples
sampleDf = sampleDataframe()
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
    print("Passed")
else:
    with open("sampleDF.csv", "r") as sampleFile:
        with open("xenaDF.csv", "r") as xenaDfFile:
            # if they are not equal then output diff of both files
            sys.stdout.writelines(difflib.unified_diff(sampleFile.readlines(), xenaDfFile.readlines(),
                                                       fromfile="sampleDF.csv", tofile="xenaDF.csv"))


