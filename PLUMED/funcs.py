import os
import urllib.request
import json
import nglview as nv
import MDAnalysis as mda


def view_nucl(*args,gui=False):
    nuclMD=mda.Universe(*args)
    #prot = nuclMD.select_atoms("protein")
    show=nv.show_mdanalysis(nuclMD,gui=gui)
    show.representations = [
    {"type": "cartoon", "params": {
        "sele": ":A :E", "color": 0x020AED,"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":B :F", "color": "green","aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":C :G", "color": 0xE0F705,"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": ":D :H", "color": 0xCE0000,"aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "cartoon", "params": {
        "sele": "nucleic", "color": "grey","aspectRatio":2, "radius":1.5,"radiusSegments":1,"capped":1
    }},
    {"type": "base", "params": {
        "sele": "nucleic", "color": "grey",
    }},
    ]
    show.camera = 'orthographic'
    return show


def get_files_from_git(gitapiurl,savefoldername):
    os.mkdir(savefoldername)
    json_url = urllib.request.urlopen(gitapiurl)
    data = json.loads(json_url.read())
    for d in data:
        if(d['type']=='file'):
            print("Downloading "+os.path.join(savefoldername,d['name']))
            urllib.request.urlretrieve(d['download_url'],os.path.join(savefoldername,d['name']))
        if(d['type']=='dir'):
            get_files_from_git(d['url'],os.path.join(savefoldername,d['name']))

  
