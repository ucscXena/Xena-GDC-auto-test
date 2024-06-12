import os
import sys
import logging
import requests
import json
import subprocess
import pandas
import numpy
import difflib


logger = logging.getLogger(__name__)


def round_ForNans(x):
    if pandas.notna(x):
        return numpy.format_float_scientific(x, precision=10)
    else:
        return numpy.nan


def vaf(t_alt_count, t_depth):
    
    return t_alt_count / t_depth


def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))
    logger.info("Downloading from GDC: ")
    subprocess.run(
        ["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
         "Content-Type: application/json", "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"])
    os.system("mkdir -p gdcFiles")
    os.system("tar -xzf gdcFiles.tar.gz -C gdcFiles")


def getXenaSamples(xenaFile):  # get all samples from the xena matrix
    xenaMatrix = pandas.read_csv(xenaFile, sep="\t")
    sampleList = list(xenaMatrix["sample"].unique())
    
    return sampleList


def getAllSamples(projectName, experimentalStrategy):
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
                        "Aliquot Ensemble Somatic Variant Merging and Masking"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        "simple nucleotide variation"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [
                        "Masked Somatic Mutation"
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
            },
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": [
                        "open"
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


def dataTypeSamples(projectName, experimentalStrategy, samples):
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
                        "Aliquot Ensemble Somatic Variant Merging and Masking"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": "simple nucleotide variation"
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": [
                        "Masked Somatic Mutation"
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
    xenaDF["dna_vaf"] = xenaDF["dna_vaf"].apply(round_ForNans)
    
    return xenaDF


def nonEmptySamples(sampleDict):
    nonEmpty = []
    allSampleNames = []
    for sample in sampleDict:
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        normalSampleName = sample[:sample.index(".")]
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)
        if (len(sampleDataDF.index) == 0):
            continue
        if (normalSampleName not in nonEmpty):
            nonEmpty.append(normalSampleName)
            allSampleNames.append(sample)

    return nonEmpty, allSampleNames


def sampleDataframe(sampleNames, sampleDict):
    dataFrame = pandas.DataFrame()
    # Create a dataframe for all the samples retrieved
    for sample in sampleNames:
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        normalSampleName = sample[:sample.index(".")]
        # Create data frame for sample data
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)
        sampleDataDF.rename(columns={'Hugo_Symbol': 'gene'}, inplace=True)
        sampleDataDF.rename(columns={'Chromosome': 'chrom'}, inplace=True)
        sampleDataDF.rename(columns={'Start_Position': 'start'}, inplace=True)
        sampleDataDF.rename(columns={'End_Position': 'end'}, inplace=True)
        sampleDataDF.rename(columns={'Reference_Allele': 'ref'}, inplace=True)
        sampleDataDF.rename(columns={'Tumor_Seq_Allele2': 'alt'}, inplace=True)
        sampleDataDF.rename(columns={'HGVSp_Short': 'Amino_Acid_Change'}, inplace=True)
        sampleDataDF.rename(columns={'Consequence': 'effect'}, inplace=True)
        sampleDataDF["dna_vaf"] = sampleDataDF.apply(lambda x: vaf(x.t_alt_count, x.t_depth), axis=1)
        sampleDataDF = sampleDataDF.loc[:,
                       ["gene", "chrom", "start", "end", "ref", "alt", "Tumor_Sample_Barcode", "Amino_Acid_Change",
                        "effect",
                        "callers", "dna_vaf"]]
        sampleDataDF = sampleDataDF[
            ["gene", "chrom", "start", "end", "ref", "alt", "Tumor_Sample_Barcode", "Amino_Acid_Change", "effect",
             "callers", "dna_vaf"]]
        sampleDataDF.insert(0, "sample", normalSampleName)
        dataFrame = pandas.concat([dataFrame, sampleDataDF])
    dataFrame["dna_vaf"] = dataFrame["dna_vaf"].apply(round_ForNans)
    
    return dataFrame


def main(projectName, xenaFilePath, dataType):
    experimentalStrategyDict = {
        "somaticmutation_wxs": "WXS",
        "somaticmutation_targeted": "Targeted Sequencing"
    }
    experimentalStrategy = experimentalStrategyDict[dataType]
    xenaSamples = getXenaSamples(xenaFilePath)
    allSamples = getAllSamples(projectName, experimentalStrategy)
    sampleDict, seenSamples = dataTypeSamples(projectName, experimentalStrategy, allSamples)
    xenaDF = xenaDataframe(xenaFilePath)
    fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
    downloadFiles(fileIDS)
    validSamples, sampleNames = nonEmptySamples(sampleDict)
    validSamples = list(set(validSamples))
    if sorted(validSamples) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(validSamples)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in validSamples if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in validSamples]}")
        exit(1)
    # create dataframe for samples
    sampleDf = sampleDataframe(sampleNames, sampleDict)
    xenaDF.sort_values(by=sorted(xenaDF), inplace=True)
    # sort sample dataframe as well
    sampleDf.sort_values(by=sorted(xenaDF), inplace=True)
    # then reset index ordering for each one
    xenaDF.reset_index(inplace=True, drop=True)
    sampleDf.reset_index(inplace=True, drop=True)
    with open("sampleDF.csv", "w") as sampleFile:
        sampleDf.to_csv(sampleFile)
    with open("xenaDF.csv", "w") as xenaDfFile:
        xenaDF.to_csv(xenaDfFile)
    try:
        logger.info("Testing in progress ...")
        pandas.testing.assert_frame_equal(sampleDf, xenaDF, check_dtype=False)
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return 'PASSED'
    except AssertionError:
        logger.info("[{}] test failed for [{}].".format(dataType, projectName))
        logger.info("Diff file is being generated with unequal values.")
        with open("sampleDF.csv", "r") as sampleFile:
            with open("xenaDF.csv", "r") as xenaDfFile:
                # if they are not equal then output diff of both files
                with open("diff.txt", "w") as diffFile:
                    diffFile.writelines(difflib.unified_diff(sampleFile.readlines(), xenaDfFile.readlines(),
                                                            fromfile="sampleDF.csv", tofile="xenaDF.csv"))
        return 'FAILED'
