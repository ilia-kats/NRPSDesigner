import json
import xml.etree.ElementTree as x

from dajaxice.decorators import dajaxice_register

@dajaxice_register
def submitNrp(request, nrp):
    parsed = json.loads(nrp)
    nrpxml = x.Element('nrp')
    for monomer in parsed:
        monomerel = x.SubElement(nrpxml, 'monomer')
        monomerid = x.SubElement(monomerel, 'id')
        monomerid.text = monomer
    return json.dumps({'xml':x.tostring(nrpxml, "utf8")})
