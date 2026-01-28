
import pathlib
import sys
import bpy
import bmesh
import random


import os

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
    sys.path.append(dir )

import list_objs_names
import organ

# this next part forces a reload in case you edit the source after you first start the blender session
import imp
imp.reload(organ)
imp.reload(list_objs_names)

import functions as fn

import imp
imp.reload(fn)


#original_mesh_obj = bpy.context.active_object
#nameObj = "12200_Skin_surface"
##nameObj = "2300_Wrists_and_hand_bones_spongiosa"
#holes = getHoles(nameObj)

#fprompt.close()

#fprompt = open("bdr_prompt.txt","w")
#fprompt.write(str(pathlib.Path().resolve())) 
#fprompt.write('\n')
#from scipy.spatial import Delaunay
#obj = bpy.data.objects['Cube'] #Using default cube as example
#me = bpy.context.object.data

scene = bpy.context.scene


#objects = bpy.data.objects
#for obj in objects:
#    print (obj.name,obj.location)

#list_objs = [obj for obj in scene.objects if obj.name.startswith("Skin")]

#ff = open("ff.txt","w")

#list_objs = [obj for obj in scene.objects]

#for obj in list_objs:
#    ff.write('"'+obj.name+'": ["'+  obj.name+ '"]\n')
#
#ff.close()
    
#    obj.name = obj.data.materials.items()[0][0].split(".")[0]

#original_mesh_obj = bpy.context.active_object
#nameObj = "2300_Wrists_and_hand_bones_spongiosa"
#nameObj = "12200_Skin_surface"
#holes = getHoles(nameObj)


#holes = getHoles("300_ET1_8um")
#fprompt.write("Holes of Et1\n")
#for h in holes:
 # fprompt.write(str(h)+"\n")
#fprompt.close()


list_objs_names = list_objs_names.list_objs_names
#listob = list_objs_nameps[0]

polydir = 'C:\sblunier3\Share\poly/'

for nameobj, listorg in list_objs_names.items():
 
    bpy.ops.object.select_all(action='SELECT')
 
    splitted = []
    org = organ.organ(listorg["out"])
    if listorg["in"] != None:
        for k,o in listorg["in"].items(): 
            list_o = [k+ol for ol in o]
            #if org.outermesh == "10600_Muscle": splitted.append(fn.split_organ(k))  
            org.addinner(k,list_o)
 
    obs = []
    all_nodes = {}
    all_tet = []

    inode = 1
    avpt = [0,0,0]
    
    listob = [bpy.context.scene.objects[org.outermesh]]
    #listob = [org.outermesh]
    
    #fprompt = open("promptd.txt","a")
    for mesh in org.splt_innermeshes:
    #for mesh in org.innermeshes:
        #fprompt.write(mesh.name + "\n")
        listob.append(mesh)
    #fprompt.close()
         
    #fprompt = open("promptb.txt","w")
    #for ob in listob:
    #    fprompt.write(ob.name+"\n")
    #fprompt.close()    
    
    #for l in listob: fprompt.write(l+" ")
    for obj in listob:
        #obs.append(scene.objects[ob])
              
        #scene.objects[obj].hide_set(False)
        obj.hide_set(False)
        
        bm = bmesh.new()
        #fprompt.write(obj+"\n")
        #bm.from_mesh(scene.objects[obj].data)   
        bm.from_mesh(obj.data)   
        bm.faces.ensure_lookup_table()
        
        maxvert = 0
        for v in bm.verts:
            all_nodes[v.index+inode] = [v.co.x,v.co.y,v.co.z]
            #fprompt.write(str(v.index+inode)+" ")
            #fprompt.write(str(v.co.x)+" ")
                #fprompt.write(str(v.co.y)+" ")
                #fprompt.write(str(v.co.z)+"\n")
            if maxvert<v.index: maxvert = v.index
            
        if inode==1:
            for k,node in all_nodes.items():
                for i in range(len(node)): avpt[i] = node[i]
            for i in range(len(avpt)): avpt[i] = avpt[i]/len(all_nodes)
        
        #fprompt.write("\n")
        for f in bm.faces:
            all_tet.append([])
            for v in f.verts: 
                all_tet[-1].append(v.index+inode)
                #fprompt.write(str(v.index+inode)+" ")
            #fprompt.write("\n") 
        
        inode = inode+maxvert+1
        
        #if org.outermesh == "10600_Muscle": 
        #    if len(splitted)>=2:
    for s in splitted:
        fn.join_objects(s)

    fpoly = open(polydir+nameobj+".poly","w")
    fpoly.write(str(len(all_nodes))+" 3\n")
    for k,v in all_nodes.items():
        fpoly.write(str(k-1))
        for val in v:
           fpoly.write(" "+str(val))
        fpoly.write("\n")
    fpoly.write("\n")

    # Add tets of organ
    fpoly.write(str(len(all_tet))+" 0 0\n")
    for tet in all_tet:
        fpoly.write("1\n")
        fpoly.write(str(len(tet))+" ")
        for t in tet:
            fpoly.write(str(t-1)+" ")
        fpoly.write("\n")
    fpoly.write("\n")

    # Add hole in organ
    
    holes = org.holes
    
    #fprompt = open("prompt.txt","a")
    #fprompt.write(nameobj+"\n")
    #for h in holes: fprompt.write(str(h)+" jj\n")
    #fprompt.close()
    
    fpoly.write(str(len(holes))+"\n")
    for ih,h in enumerate(holes):
        fpoly.write(str(ih+1)+"x ")
        #fprompt.write(str(h)+"\n")
        #fprompt.close()
        for c in h: fpoly.write(str(c)+" ")
        fpoly.write("\n")

    fpoly.write("\n\n0\n")
    fpoly.close()
    
    org.join_components()
    bpy.ops.object.select_all(action='SELECT')


#fprompt.close()