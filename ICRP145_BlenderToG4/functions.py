import mathutils
import bpy
import bmesh
import random
import numpy as np

def is_point_inside_mesh(mesh_obj, point, num_rays=10):

    directions = [mathutils.Vector((random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1))).normalized() for _ in range(num_rays)]

    for direction in directions:
        result, pt1, normal, index = mesh_obj.ray_cast(point,direction)
        resultneg, pt2, normalneg, indexneg = mesh_obj.ray_cast(point,-direction)
        if not (result and resultneg): return False   
        
    return True

def is_object_inside_mesh(mesh_obj, obj):
    
    bm = bmesh.new()
    fprompt = open("prompt.txt",'a')
    fprompt.write(obj.name+"\n")
    fprompt.close()
    bm.from_mesh(obj.data)   
    bm.faces.ensure_lookup_table()

    
    maxvert = 0
    for v in bm.verts:
        if is_point_inside_mesh(mesh_obj,mathutils.Vector([v.co.x,v.co.y,v.co.z])): return True
        else: return False
     

def getIntersections(mesh_obj,point,direction,tol=0.001):
    list_intersection = []
    #bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)
    while True:
        result, point, normal, index = mesh_obj.ray_cast(point+direction*tol,direction)
        if result: 
            list_intersection.append(point)
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)
        else: break
        #_intersection = n_intersection+1
    return list_intersection

def is_point_inside_mesh2(mesh_obj, point, num_rays=10):

    directions = [mathutils.Vector((random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1))).normalized() for _ in range(num_rays)]

    for direction in directions:
        if len(getIntersections(mesh_obj,point,direction))%2==0: return False
        if len(getIntersections(mesh_obj,point,-direction))%2==0: return False

    return True


def find_internal_point(mesh_obj,num_rays=10):
    
    #scene = bpy.context.scene

    #depsgraph = bpy.context.evaluated_depsgraph_get()
    #
    #obj_eval = mesh_obj.evaluated_get(depsgraph)
    
    
    # Ensure the object is a mesh
    if mesh_obj.type != 'MESH':
        raise TypeError("The object must be a mesh")
    
    #fprompt.write("obj_eval = "+str(obj_eval))
    directions = [mathutils.Vector((random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1))).normalized() for _ in range(num_rays)]
   
    i=0
    
    bm = bmesh.new()
    bm.from_mesh(mesh_obj.data)   
    bm.faces.ensure_lookup_table()
    
    fprompt = open("prompt.txt",'a')
    fprompt.write(mesh_obj.name+"\n")
    fprompt.close()
    #bm.from_mesh(obj.data)   
    #bm.faces.ensure_lookup_table()

    
    n = len(bm.verts)  
    while i<1000:
        i1 = np.random.randint(n)
        i2 = np.random.randint(n)
        while i1 == i2: 
            i2 = np.random.randint(n)        
            
        for iv, v in enumerate(bm.verts):
            if iv == i1: 
                pt1 = mathutils.Vector([v.co.x,v.co.y,v.co.z])
                continue
            if iv == i2: 
                pt2 = mathutils.Vector([v.co.x,v.co.y,v.co.z])
                continue
                    
        haspoint = False

        maxdistance = 0


        #fprompt = open("testHoles.txt","w")
        for direction in directions:
            # Cast the ray
            point = (pt1 + pt2)/2.0
            #result, pt1, normal, index, hit_obj, matrix = scene.ray_cast(depsgraph, point, direction)
            #resultneg, pt2, normalneg, indexneg, hit_objneg, matrixneg = scene.ray_cast(depsgraph, point, -direction)
            #result, pt1, normal, index = mesh_obj.ray_cast(point,direction)
            
            result, pt1, normal1, index1 = mesh_obj.ray_cast(point,direction)
            #fprompt.write(str(result1)) 
            #fprompt.write("\n")
            
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=pt1)
            
            #getIntersections(mesh_obj,pt1,direction)
            
            #break
            resultneg, pt2, normalneg, indexneg = mesh_obj.ray_cast(point,-direction)
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=pt1+direction*0.01)
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=pt11)
            
            #break
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)
            # Make sure the selected point is in the meshed structure
            if not (result and resultneg): 
                haspoint = False
                break
            else:
                haspoint = True
            #point = (locationneg+location)/2.0
        
        #fprompt.write(str(len(getIntersections(mesh_obj,point,direction))))
        
        #fprompt.close()
        #if not is_point_inside_mesh(mesh_obj,point):
        #    haspoint = False 
    
        
    
        for direction in directions:
            if len(getIntersections(mesh_obj,point,direction))%2==0:
                haspoint = False
                break
        i = i+1 
        if haspoint: 
            #bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)
            return point
    
        #if is_point_inside_mesh(mesh_obj,point): 
        #    bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)
        #    return point
    
          
     
    raise ValueError("No point found in object "+ str(mesh_obj.name))    
#'''

def split_organ(nameObjout):

    #nameObjout = "10100_Lymphatic_nodes_thoracic"
    original_mesh_objout = bpy.context.scene.objects[nameObjout]
    bpy.data.objects[nameObjout].select_set(True)
    bpy.context.view_layer.objects.active = bpy.context.scene.objects[nameObjout]

    # Separate the connected components
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.separate(type='LOOSE')
    bpy.ops.object.mode_set(mode='OBJECT')
    
    bpy.ops.object.select_all(action='DESELECT')
    
    return bpy.context.selected_objects

def join_objects(objects):
    """
    Join a list of objects into one and name the resulting object after the first object in the list.
    
    Parameters:
    objects (list): A list of bpy.types.Object to join.
    """
    #if not objects:
     #   raise ValueError("The objects list is empty")
        
    if len(objects)<=1: return

    fprompt = open("split_objects.txt","a")
    for obj in objects:
        fprompt.write(obj.name+"\n")
    fprompt.close()


    # Deselect all objects
    bpy.ops.object.select_all(action='DESELECT')
    
    # Select the objects to join
    for obj in objects:
        obj.select_set(True)
    
    # Make the first object active
    bpy.context.view_layer.objects.active = objects[0]
    
    # Join the objects
    bpy.ops.object.join()
    
    # Rename the resulting object
    joined_object = bpy.context.view_layer.objects.active
    joined_object.name = objects[0].name
    joined_object.data.name = objects[0].name
    joined_object.active_material.name = objects[0].name


def getHoles(nameObj):
    
    #fprompt = open("prompt2.txt","a")
    #fprompt.write(nameObj)
    
    bpy.ops.object.select_all(action='DESELECT')
    #bpy.data.objects["Cube"].select_set(True)
    bpy.context.view_layer.objects.active = bpy.data.objects[nameObj]
    original_mesh_obj = bpy.context.scene.objects[nameObj]
    bpy.data.objects[nameObj].select_set(True)
    #print (nameObj)
    bpy.context.view_layer.objects.active = bpy.context.scene.objects[nameObj]

    # Separate the connected components
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.separate(type='LOOSE')
    bpy.ops.object.mode_set(mode='OBJECT')
        
    connected_components = bpy.context.selected_objects
    #fprompt.write(str(len(connected_components))+"\n")
    #for c in connected_components:[
     #   fprompt.write(c.name+"\n")
    #fprompt.write("After connected components\n")
   
    ## Find an internal point for each connected component
    internal_points = []
    i=0
    #split_objects = []
    for component in connected_components:
        point = find_internal_point(component) 
        internal_points.append(point)
        #break
        #split_objects.append(component.copy())
        #fprompt.write(str(point)+"\n")
        #bpy.ops.object.empty_add(type='PLAIN_AXES', location=point)

    #if len(connected_components)>0: join_objects(connected_components)

    #fprompt = open("split_objects.txt","a")
    #fprompt.write("in Get Holes")
    #for obj in connected_components:
    #    fprompt.write(obj.name+"\n")
    #fprompt.close()
 
    bpy.ops.object.select_all(action='DESELECT')
    
    #fprompt = open("prompta.txt","a")
    #for spltobj in split_objects: fprompt.write(spltobj.name)
    #fprompt.close()
    return internal_points, connected_components

#getHoles("10300_Lymphatic_nodes_trunk")