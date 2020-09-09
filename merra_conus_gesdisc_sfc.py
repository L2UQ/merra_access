# Query and download MERRA 2 meteorological fields from GES-DISC
# Radiation diagnostics: TS, CLDTOT

# STEP 1
import sys
import json
import urllib3
import certifi
import requests
from time import sleep
from http.cookiejar import CookieJar
import urllib.request
from urllib.parse import urlencode

# STEP 2
# Create a urllib PoolManager instance to make requests.
http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',ca_certs=certifi.where())
# Set the URL for the GES DISC subset service endpoint
url = 'https://disc.gsfc.nasa.gov/service/subset/jsonwsp'

# STEP 3
# This method POSTs formatted JSON WSP requests to the GES DISC endpoint URL
# It is created for convenience since this task will be repeated more than once
def get_http_data(request):
    hdrs = {'Content-Type': 'application/json',
            'Accept'      : 'application/json'}
    data = json.dumps(request)       
    r = http.request('POST', url, body=data, headers=hdrs)
    response = json.loads(r.data.decode('utf-8'))   
    # Check for errors
    if response['type'] == 'jsonwsp/fault' :
        print('API Error: faulty %s request' % response['methodname'])
        sys.exit(1)
    return response

# STEP 4
# Define the parameters for the data subset
# Time range might need to extend 1 day beyond desired
product = 'M2I1NXASM_V5.12.4' 
varNames =['TS','PS']
minlon = -140.0
maxlon = -40.0
minlat = 20.0
maxlat = 60.0
begTime = '2019-04-30'
endTime = '2019-05-01'
begHour = '09:00'
endHour = '09:00'
doMean = False

# Subset only the mandatory pressure levels (units are hPa)
# 1000 925 850 700 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1 
#dimName = 'lev'
#dimVals = [1,4,7,13,17,19,21,22,23,24,25,26,27,29,30,31,32,33,35,36,37]
# Construct the list of dimension name:value pairs to specify the desired subset
#dimSlice = []
#for i in range(len(dimVals)):
#    dimSlice.append({'dimensionId': dimName, 'dimensionValue': dimVals[i]})

# STEP 5
# Construct JSON WSP request for API method: subset
subset_request = {
    'methodname': 'subset',
    'type': 'jsonwsp/request',
    'version': '1.0',
    'args': {
        'role'  : 'subset',
        'box'   : [minlon, minlat, maxlon, maxlat],
        'start' : begTime,
        'end'   : endTime,
        'diurnalFrom': begHour, 
        'diurnalTo': endHour,
        'diurnalMean': doMean,
        'data': [{'datasetId': product,
                  'variable' : varNames[0]
                 },
                 {'datasetId': product,
                  'variable' : varNames[1]
                 }]
    }
}

# STEP 6
# Submit the subset request to the GES DISC Server
response = get_http_data(subset_request)
# Report the JobID and initial status
myJobId = response['result']['jobId']
print('Job ID: '+myJobId)
print('Job status: '+response['result']['Status'])

# STEP 7
# Construct JSON WSP request for API method: GetStatus
status_request = {
    'methodname': 'GetStatus',
    'version': '1.0',
    'type': 'jsonwsp/request',
    'args': {'jobId': myJobId}
}

# Check on the job status after a brief nap
while response['result']['Status'] in ['Accepted', 'Running']:
    sleep(5)
    response = get_http_data(status_request)
    status  = response['result']['Status']
    percent = response['result']['PercentCompleted']
    print ('Job status: %s (%d%c complete)' % (status,percent,'%'))
if response['result']['Status'] == 'Succeeded' :
    print ('Job Finished:  %s' % response['result']['message'])
else : 
    print('Job Failed: %s' % response['fault']['code'])
    sys.exit(1)

# STEP 8 
# Construct JSON WSP request for API method: GetResult
batchsize = 20
results_request = {
    'methodname': 'GetResult',
    'version': '1.0',
    'type': 'jsonwsp/request',
    'args': {
        'jobId': myJobId,
        'count': batchsize,
        'startIndex': 0
    }
}

# Retrieve the results in JSON in multiple batches 
# Initialize variables, then submit the first GetResults request
# Add the results from this batch to the list and increment the count
results = []
count = 0 
response = get_http_data(results_request) 
count = count + response['result']['itemsPerPage']
results.extend(response['result']['items']) 

# Increment the startIndex and keep asking for more results until we have them all
total = response['result']['totalResults']
while count < total :
    results_request['args']['startIndex'] += batchsize 
    response = get_http_data(results_request) 
    count = count + response['result']['itemsPerPage']
    results.extend(response['result']['items'])
        
# Check on the bookkeeping
print('Retrieved %d out of %d expected items' % (len(results), total))

# STEP 9 
# Sort the results into documents and URLs
docs = []
urls = []
for item in results :
    try:
        if item['start'] and item['end'] : urls.append(item) 
    except:
        docs.append(item)
# Print out the documentation links, but do not download them
print('\nDocumentation:')
for item in docs : print(item['label']+': '+item['link'])

print(urls)
# STEP 10 
# Use the requests library to submit the HTTP_Services URLs and write out the results.
print('\nHTTP_services output:')
for item in urls :
    URL = item['link'] 
    result = requests.get(URL)
    try:
        result.raise_for_status()
        outfn = item['label']
        f = open(outfn,'wb')
        f.write(result.content)
        f.close()
        print(outfn)
    except:
        print('Error! Status code is %d for this URL:\n%s' % (result.status.code,URL))
        print('Help for downloading data is at https://disc.gsfc.nasa.gov/data-access')
print('OK')

