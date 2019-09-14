# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:04:19 2018

@author: Phil

This code is copied from the website
http://idl2python.blogspot.com/2010/10/save-and-restore-2.html


"""
#
#
#!!!!!!!!!!!!!!!!!!!!!
#this is writen for python 2.7, needs updating for python 3
#
#
#
#
#def save(file,**kwargs):
#    """
#    Save the value of some data in a file.
#    Usage: save('misdatos.pypic',a=a,b=b,test=test)
#    """
##   import cPickle
#    f=open(file,"wb")
#    pickle.dump(kwargs,f,protocol=2)
#    f.close
#
#def restore(file):
#    """
#    Read data saved with save function.
#    Usage: datos = restore('misdatos.pypic')
#    """
##   import cPickle
#    f=open(file,"rb")
#    result = pickle.load(f)
#    f.close
#    return result
#
