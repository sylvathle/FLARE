import pathlib
import sys
import bpy


import numpy as np

import os

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import functions as fn

import imp
imp.reload(fn)

class organ:
    outermesh = ""
    innermeshes = []
    holes = []
    
    def __init__(self,outermesh =[]):
        self.outermesh = outermesh
        self.holes = []
        self.innermeshes = []
        self.list_components = []
        self.splt_innermeshes = []
        #self.splt_innerexcluded = []
        #self.innermeshes = innermeshes
        #for inmesh in innermeshes:
        #    self.innermeshes.append(inmesh)
            #fprompt.write(inmesh+"\n")
            #fprompt.close()
        #    hol = getHoles(inmesh)
        #    for h in hol: self.holes.append(h)
    
    def addinner(self,inner,excluded=[]):
        self.innermeshes.append(inner)
        #fprompt = open("promptc.txt","a")
        
        hol, splt_objs = fn.getHoles(inner)
        self.list_components.append(splt_objs)
        #self.splt_innerexcluded.append(exclude)
        #self.holes.append(hol)
        #hol = fn.getHoles(inner)
        #for obj in splt_objs: 
        #    if fn.is_object_inside_mesh(bpy.context.scene.objects[self.outermesh],obj):
        #       self.splt_innermeshes.append(obj)
        #self.splt_innermeshes.append(splt_objs)
        for ii,h in enumerate(hol): 
            #fprompt.write(self.outermesh+"\t")
            #fprompt.write(inner+"\t")
            #fprompt.write(splt_objs[ii].name+"\t")
            #fprompt.write(str(h)+"\t")1
            isInouterMesh = fn.is_point_inside_mesh(bpy.context.scene.objects[self.outermesh],h)
            notInMuscle = (self.outermesh!="11600_RST") and (not fn.is_point_inside_mesh2(bpy.context.scene.objects["10600_Muscle"],h)) or (inner=="10600_Muscle") 
            notInMuscle = True
            if isInouterMesh and notInMuscle:
                self.holes.append(h)
                ##if splt_objs[ii].name not in excluded:
                if splt_objs[ii].name not in excluded:      
                    self.splt_innermeshes.append(splt_objs[ii])
                    #bpy.ops.object.empty_add(type='PLAIN_AXES', location=h)
                #else: bpy.ops.object.empty_add(type='PLAIN_AXES', location=h)
                #fprompt.write("IN")
            #else: 
                #fprompt.write("OUT")
                #bpy.ops.object.empty_add(type='PLAIN_AXES', location=h)
            #fprompt.write("\n")
        #fprompt.close()
        
    def join_components(self):
        
        #fprompt = open("split_objects.txt","a")
    
        #for obj in objects:
        #fprompt.write("In join components \n")
        #fprompt.close()
    
        for objects in self.list_components:
            if len(objects)>1: fn.join_objects(objects)
