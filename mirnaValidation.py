import requests
import logging
import json
import tarfile
import subprocess
import pandas
import numpy
import sys
import os
from math import log10, floor


logger = logging.getLogger(__name__)


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
    logger.info("Downloading from GDC: ")
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
            if sampleName in mirnaSamplesDict:
                mirnaSamplesDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]
            else:
                mirnaSamplesDict[sampleName] = {}
                mirnaSamplesDict[sampleName][caseDict["file_id"]] = caseDict["file_name"]
    
    return mirnaSamplesDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t", index_col=0)
    for column in xenaDF:
        xenaDF[column] = xenaDF[column].apply(round_ForNans)
    return xenaDF

def mirnaDataframe(mirnaSamplesDict, dataType):
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
            fileName = mirnaSamplesDict[sample][fileID]
            sampleFile = "gdcFiles/" + fileID + "/" + fileName
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
    fileIDS = [fileID for sample in mirnaSamplesDict for fileID in mirnaSamplesDict[sample]]
    downloadFiles(fileIDS)
    gdcDF = mirnaDataframe(mirnaSamplesDict, dataType)
    result = compare(logger, gdcDF, xenaDF)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return 'PASSED'
    else:
        logger.info("[{}] test failed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return 'FAILED'
