#!/usr/bin/env python
from __future__ import print_function
import os, sys, glob, shutil, multiprocessing
import math


def run(filename):
    to_glob = filename[:filename.rfind('_')] + '_*.root'
    # print(to_glob)
    list_of_files = glob.glob(to_glob)
    list_of_files.sort(key=lambda x: int(x[x.rfind('_')+1:x.rfind('.')]))
    # print(list_of_files)

    inputs = []
    for f in list_of_files:
        cmd = "root -b -q -l 'run.C(\"%s\")'" % (f,)
        inputs.append(cmd)

    nCores = multiprocessing.cpu_count()
    print('total cpu cores: ', nCores)
    pool = multiprocessing.Pool(nCores*5)
    for cmd in inputs:
        pool.apply_async(os.system, args=(cmd,))
    pool.close()
    pool.join()

def organize():
    plotDir = 'tmp_plots'
    os.chmod(plotDir, 0775)
    umask_changed = False
    original_umask = os.umask(002)
    if (original_umask != 002):
        umask_changed = True

    os.chdir(plotDir)
    files = os.listdir('.')
    # print(files)
    for f in files:
        fullname, ext = os.path.splitext(f)
        if ext in ('.png', '.jpg'):
            # print(name, ext)
            ind = fullname.rfind('_')
            name = fullname[:ind]
            eventId = fullname[ind+1:]
            run, subrun, event = eventId.split('-')
            thousands = '%06d' % (math.floor(int(run) / 1000) * 1000, )
            dirname = '%s/%s/%s/%s' % (thousands, run, subrun, event)
            new_f = dirname+'/'+name+ext
            # print(f, new_f)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            os.rename(f, new_f)
            os.chmod(new_f, 0664)
    os.chdir('..')
    if (umask_changed):
        os.umask(original_umask)
    print('re-orginazing to folders ... done.')

def usage():
    print("""
    python makeplots.py [filename]
    """)

if __name__ == "__main__":
    if (len(sys.argv)<=1):
        usage()
    else:
        args = sys.argv
        args.pop(0)
        filename = args.pop(0)
        # print(filename)
        run(filename)
        organize()

