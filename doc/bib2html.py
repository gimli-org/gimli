import json
import os

from bibtexparser import load
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode


def parse_bib(fname):
    """ Read bibtex file and sort by year. """

    with open(fname) as bibfile:
        parser = BibTexParser()
        parser.customization = convert_to_unicode
        bp = load(bibfile, parser=parser)
        references = bp.get_entry_list()

    references.sort(key=lambda x: x['year'], reverse=True)

    return references


def write_html():
    db = parse_bib("gimliuses.bib")
    for entry in db:
        entry["author"] = entry["author"].replace(" and ", ", ")
        if not "journal" in entry:
            entry["journal"] = entry.pop("booktitle")
        if not "doi" in entry:
            string = ""
        else:
            doi = entry["doi"]
            string = "<a target='_blank' href='https://doi.org/%s'>%s</a>" % (
                doi, doi)
        entry["doi"] = string

    return json.dumps(db, sort_keys=True, indent=4)
