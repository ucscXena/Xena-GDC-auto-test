import json
import pandas
import requests
import logging

logger = logging.getLogger(__name__)


def getTimeAndPatientData(project):
    survivalEndpoint = "https://api.gdc.cancer.gov/analysis/survival"

    filter = {
        "op": "in",
        "content": {
            "field": "cases.project.project_id",
            "value": [project]
        }
    }
    params = {
        "filters": filter,
        "fields": "",
        "format": "json",
        "size": 200000
    }

    response = requests.post(survivalEndpoint, headers={"Content-Type": "application/json"}, json=params)
    responseJson = response.json()
    timeData = responseJson["results"][0]["donors"]
    timeDict = {}
    for entry in timeData:
        patient = entry["submitter_id"]
        osTime = entry["time"]
        timeDict[patient] = {"OS.time": osTime}
    return timeDict


def getOS(submitterIDs, survivalData):
    casesEndpoint = "https://api.gdc.cancer.gov/cases"
    filter = {
        "op": "in",
        "content": {
            "field": "submitter_id",
            "value": submitterIDs
        }
    }
    params = {
        "filters": filter,
        "fields": "demographic.vital_status,submitter_id,submitter_sample_ids",
        "format": "json",
        "size": 2000000
    }

    response = requests.post(casesEndpoint, headers={"Content-Type": "application/json"}, json=params)
    responseJson = response.json()
    statusData = responseJson["data"]["hits"]
    for entry in statusData:
        submitterID = entry["submitter_id"]
        status = 1 if entry["demographic"]["vital_status"] == "Dead" else 0
        samples = entry["submitter_sample_ids"]
        survivalData[submitterID]["OS"] = status
        survivalData[submitterID]["samples"] = samples
    return survivalData


def getKeepSamples(project):
    filesEndpoint = "https://api.gdc.cancer.gov/files/"
    fileFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        project
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "access",
                    "value": [
                        'open'
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        'transcriptome profiling', 'proteome profiling', 'dna methylation', 'copy number variation',
                        'simple nucleotide variation'
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(fileFilter),
        "fields": 'data_category,cases.samples.submitter_id,cases.samples.tissue_type',
        "format": "json",
        "size": 2000000
    }

    response = requests.post(filesEndpoint, json=params,
                             headers={"Content-Type": "application/json"})
    response = response.json()["data"]["hits"]

    keepSamples = []
    for file in response:
        samples = file["cases"][0]["samples"]
        # Get samples associated with transcriptome, proteome, and/or methylation file(s)
        if file['data_category'] in ['Transcriptome Profiling', 'Proteome Profiling', 'DNA Methylation']:
            for sample in samples:
                submitter_id = sample['submitter_id']
                if submitter_id not in keepSamples:
                    keepSamples.append(submitter_id)
        # Get 'Tumor' samples associated with copy number variation and/or simple nucleotide variation file(s)
        if file['data_category'] in ['Copy Number Variation', 'Simple Nucleotide Variation']:
            for sample in samples:
                submitter_id = sample['submitter_id']
                if sample['tissue_type'] == 'Tumor' and submitter_id not in keepSamples:
                    keepSamples.append(submitter_id)
    return keepSamples


def processSamples(survivalData, keepSamples):
    for patient in survivalData:
        survivalData[patient]["samples"] = [sample for sample in survivalData[patient]["samples"] if sample in
                                            keepSamples]
    return survivalData


def formatData(survivalData):
    formattedData = [
        {
            "sample": sample,
            "OS.time": survivalData[patient]["OS.time"],
            "OS": survivalData[patient]["OS"],
            "_PATIENT": patient
         }
        for patient in survivalData for sample in survivalData[patient]["samples"]]
    return formattedData


def compare(logger, gdcDF, xenaDF):
    failed = []
    sampleNum = 1
    total = len(gdcDF.index)
    for sample in gdcDF["sample"]:
        xenaRow = xenaDF.loc[sample]
        gdcRow = gdcDF.loc[sample]
        if not xenaRow.equals(gdcRow):
            status = "[{:d}/{:d}] Sample: {} - Failed"
            logger.info(status.format(sampleNum, total, sample))
            failed.append('{} ({})'.format(sample, sampleNum))
        else:
            status = "[{:d}/{:d}] Sample: {} - Passed"
            logger.info(status.format(sampleNum, total, sample))
        sampleNum += 1

            # print("--------------")
            # print("Error: Xena Row and GDC Row are not equal\n")
            # print("Xena Row")
            # print(xenaRow, end="\n\n")
            # print("GDC Row")
            # print(gdcRow)
            # print("--------------\n\n")
    return failed


def main(projectName, xenaFilePath, dataType):
    timeData = getTimeAndPatientData(projectName)
    submitterIDs = [submitterId for submitterId in timeData]
    survivalData = getOS(submitterIDs, timeData)
    keepSamples = getKeepSamples(projectName)
    survivalData = processSamples(survivalData, keepSamples)
    formattedData = formatData(survivalData)
    gdcDF = pandas.DataFrame(formattedData)
    xenaDF = pandas.read_csv(xenaFilePath, sep="\t")

    gdcDF.sort_values(by=["sample"], inplace=True)
    xenaDF.sort_values(by=["sample"], inplace=True)
    gdcDF.reset_index(inplace=True, drop=True)
    xenaDF.reset_index(inplace=True, drop=True)
    gdcDF.set_index("sample", inplace=True, drop=False)
    xenaDF.set_index("sample", inplace=True, drop=False)

    xenaSamples = xenaDF["sample"].tolist()
    gdcSamples = gdcDF["sample"].tolist()
    if not gdcSamples == xenaSamples:
        logger.info("ERROR: Samples retrieved from the GDC do not match those found in Xena matrix.")
        logger.info(f"Number of samples from the GDC: {len(gdcSamples)}")
        logger.info(f"Number of samples in Xena matrix: {len(xenaSamples)}")
        logger.info(f"Samples from GDC and not in Xena: {[x for x in gdcSamples if x not in xenaSamples]}")
        logger.info(f"Samples from Xena and not in GDC: {[x for x in xenaSamples if x not in gdcSamples]}")
        exit(1)

    logger.info("Testing in progress ...")

    result = compare(logger, gdcDF, xenaDF)
    if len(result) == 0:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        return "PASSED"
    else:
        logger.info("[{}] test passed for [{}].".format(dataType, projectName))
        logger.info("Samples failed: {}".format(result))
        return "FAILED"
