import os
import logging
import requests
import json
import subprocess
import pandas
import numpy
import difflib
import hashlib

logger = logging.getLogger(__name__)


def md5sum(filePath):
    md5_hash = hashlib.md5()

    with open(filePath, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()


def existing_md5sums(logger, projectName, dataType, sampleDict):
    existingFiles = [file for file in os.listdir(f"gdcFiles/{projectName}/{dataType}") if not file.startswith(".")]
    existingMd5sumFileDict = {md5sum(f"gdcFiles/{projectName}/{dataType}/{file}"): file for file in existingFiles}
    gdcMd5sumFileDict = {sampleDict[sample][fileID]["md5sum"]: fileID for sample in sampleDict for fileID in
                         sampleDict[sample]}
    logger.info(f"{len(gdcMd5sumFileDict)} files found from the GDC for {dataType} data for {projectName}")
    logger.info(f"{len(existingMd5sumFileDict)} files found at gdcFiles/{projectName}/{dataType}")
    fileIdDict = {innerKey: value
                  for outerDict in sampleDict.values()
                  for innerKey, value in outerDict.items()}
    filesNeededToUpdate = {gdcMd5sumFileDict[md5sum]: fileIdDict[gdcMd5sumFileDict[md5sum]]["file_name"] for md5sum in
                           gdcMd5sumFileDict if md5sum not in existingMd5sumFileDict}
    logger.info(f"{len(filesNeededToUpdate)} files needed to update")
    return filesNeededToUpdate


def round_ForNans(x):
    if pandas.notna(x):
        return numpy.format_float_scientific(x, precision=8)
    else:
        return numpy.nan


def vaf(t_alt_count, t_depth):
    return t_alt_count / t_depth


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
        "fields": "cases.samples.submitter_id,cases.samples.tissue_type,file_id,file_name,md5sum",
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
                dataTypeDict[sampleName + "." + str(seenDict[sampleName])] = {
                    caseDict["file_id"]: {
                        "file_name": caseDict["file_name"],
                        "md5sum": caseDict["md5sum"]}
                }

    return dataTypeDict, list(set(seenSamples))


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t")
    xenaDF["dna_vaf"] = xenaDF["dna_vaf"].apply(round_ForNans)

    return xenaDF


def nonEmptySamples(sampleDict, projectName, dataType):
    nonEmpty = []
    allSampleNames = []
    for sample in sampleDict:
        for fileID in sampleDict[sample]:
            fileName = sampleDict[sample][fileID]["file_name"]
            sampleFile = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            normalSampleName = sample[:sample.index(".")]
            sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)
            if len(sampleDataDF.index) == 0:
                allSampleNames.append(normalSampleName)
                continue
            if sample not in nonEmpty:
                nonEmpty.append(sample)
                allSampleNames.append(normalSampleName)

    return nonEmpty, allSampleNames


def sampleDataframe(sampleDict, nonEmpty, projectName, dataType):
    dataFrame = pandas.DataFrame()
    # Create a dataframe for all the samples retrieved
    for sample in sampleDict:
        for fileID in sampleDict[sample]:
            fileName = sampleDict[sample][fileID]["file_name"]
            sampleFile = "gdcFiles/{}/{}/{}".format(projectName, dataType, fileName)
            normalSampleName = sample[:sample.index(".")]
            # Create data frame for sample data
            if sample in nonEmpty:
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
            else:
                noMutation = {
                    'gene': numpy.nan,
                    'chrom': numpy.nan,
                    'start': -1,
                    'end': -1,
                    'ref': numpy.nan,
                    'alt': numpy.nan,
                    'Amino_Acid_Change': numpy.nan,
                    'effect': numpy.nan,
                }
                sampleDataDF = pandas.DataFrame([noMutation])
                sampleDataDF.insert(0, 'sample', normalSampleName)
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
    if os.path.isdir(f"gdcFiles/{projectName}/{dataType}"):
        fileIDs = existing_md5sums(logger, projectName, dataType, sampleDict)
    else:
        fileIDs = [fileID for sample in sampleDict for fileID in sampleDict[sample]]
        logger.info(f"{len(fileIDs)} files found from the GDC for {dataType} data for {projectName}")
        logger.info(f"0 files found at gdcFiles/{projectName}/{dataType}")
        logger.info(f"{len(fileIDs)} files needed to download")
    if len(fileIDs) != 0:
        downloadFiles(fileIDs, projectName, dataType)
    nonEmpty, sampleNames = nonEmptySamples(sampleDict, projectName, dataType)
    nonEmpty = list(set(nonEmpty))
    sampleNames = list(set(sampleNames))
    if sorted(sampleNames) != sorted(xenaSamples):
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(sampleNames)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in sampleNames if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in sampleNames]}")
        exit(1)
    # create dataframe for samples
    sampleDf = sampleDataframe(sampleDict, nonEmpty, projectName, dataType)
    xenaDF.sort_values(by=sorted(xenaDF), inplace=True)
    # sort sample dataframe as well
    sampleDf.sort_values(by=sorted(xenaDF), inplace=True)
    # then reset index ordering for each one
    xenaDF.reset_index(inplace=True, drop=True)
    sampleDf.reset_index(inplace=True, drop=True)

    # drop empty mutation representations if there is existing data in dataframe
    noMutationIndices = sampleDf.index[sampleDf['start'] == -1].tolist()
    indicesToDrop = []
    for noMutationIndex in noMutationIndices:
        noMutationSample = sampleDf.iloc[noMutationIndex]["sample"]
        if not xenaDF[(xenaDF["sample"] == noMutationSample) & (xenaDF['start'] != -1)].empty:
            indicesToDrop.append(noMutationIndex)
    sampleDf.drop(indicesToDrop, inplace=True)
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
