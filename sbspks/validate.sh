#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH="${PYTHONPATH}:$DIR/../Django/DjangoNrps"
export DJANGO_SETTINGS_MODULE="DjangoNrps.settings"

./validate.py $@

