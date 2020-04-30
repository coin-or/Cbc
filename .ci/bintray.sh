#!/usr/bin/bash

while [[ ! `curl -utkralphs:$BINTRAY_API -X PUT -H "Content-Type: application/json" -d'{"list_in_downloads":true}' https://api.bintray.com/file_metadata/coin-or/download/blah.zip` =~ success ]]
do
    echo "Command failed... trying again in 10 seconds"
    #No sleep command in MSys bash
    coproc read -t 10 && wait "$!" || true
done
