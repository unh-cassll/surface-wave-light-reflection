"""Render the exported surface in Blender (Cycles) for an INTENSITY-only
cross-check of seapol's I channel.

Run from outside Python, inside Blender:
    blender -b --python validation/render_blender.py -- \
        validation/output/scene.json validation/output/surface.obj \
        validation/output/blender.exr

The surface is a reflection-only dielectric reflectance map: a sharp mirror
reflecting a uniform sky (color = SceneSpec.sky_rgb), weighted by the dielectric
Fresnel reflectance, non-reflected fraction absorbed.  So the rendered value is
R_Fresnel(slope) * sky -- the slope-modulated sky reflectance, the analog of
seapol's render_camera_image with water=None.

Caveats (intentional, documented):
    * Blender has no polarization -- this validates radiance only, never
      DoLP/AoP.
    * Blender's sky/dielectric model differs from seapol's, and a Glass BSDF in
      a uniform field would just return that field (flat).  The reflectance-map
      material above gives the slope structure; the comparison is qualitative
      (hue and wave pattern), not absolute radiance.
"""

import sys

import numpy as np

# args after "--"
argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
scene_json, obj_path, out_exr = argv[0], argv[1], argv[2]

# make scene.py importable regardless of Blender's cwd
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from scene import SceneSpec  # noqa: E402

import bpy  # noqa: E402
import mathutils  # noqa: E402

spec = SceneSpec.from_json(scene_json)
cb = spec.camera_basis()
n_water = spec.n_water if spec.n_water > 0 else 1.34


# --- clean scene ----------------------------------------------------------
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = "CYCLES"

# --- surface mesh ---------------------------------------------------------
# up_axis="Z", forward_axis="Y" = identity import: the OBJ is already in
# seapol's Z-up frame, so suppress Blender's default Y-up axis conversion
# (which would otherwise rotate the flat patch into a vertical sliver).
bpy.ops.wm.obj_import(filepath=obj_path, forward_axis="Y", up_axis="Z")
mesh_obj = bpy.context.selected_objects[0]
bpy.context.view_layer.objects.active = mesh_obj

# Ensure facet normals point up (+z).  An inverted/backfacing mesh makes the
# Fresnel node see grazing incidence everywhere -> reflectance ~1, washing out
# the slope structure.  Recalculate consistently, then flip if the mean normal
# points down.
bpy.ops.object.mode_set(mode="EDIT")
bpy.ops.mesh.select_all(action="SELECT")
bpy.ops.mesh.normals_make_consistent(inside=False)
bpy.ops.object.mode_set(mode="OBJECT")
mesh_obj.data.update()
nz = sum(p.normal[2] for p in mesh_obj.data.polygons) / max(len(mesh_obj.data.polygons), 1)
if nz < 0.0:
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.flip_normals()
    bpy.ops.object.mode_set(mode="OBJECT")
    mesh_obj.data.update()
    nz = sum(p.normal[2] for p in mesh_obj.data.polygons) / max(len(mesh_obj.data.polygons), 1)
print(f"blender mesh mean normal z = {nz:+.3f} (want > 0, surface faces up)")

# Reflection-only dielectric reflectance map: a sharp mirror reflecting the
# uniform sky, weighted by the dielectric Fresnel reflectance (Fresnel node,
# IOR = seawater), the non-reflected fraction absorbed (black).  With the mesh
# normals corrected to +z (above), the Fresnel node gives the true
# slope-dependent reflectance, so S0 = R_Fresnel(slope) * I_sky -- the analog
# of seapol's reflection-only render.  Blender is unpolarized: qualitative only.
mat = bpy.data.materials.new("seawater")
mat.use_nodes = True
nt = mat.node_tree
nt.nodes.clear()
out = nt.nodes.new("ShaderNodeOutputMaterial")
mix = nt.nodes.new("ShaderNodeMixShader")
mirror = nt.nodes.new("ShaderNodeBsdfGlossy")
mirror.inputs["Roughness"].default_value = 0.0
black = nt.nodes.new("ShaderNodeBsdfDiffuse")
black.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
fresnel = nt.nodes.new("ShaderNodeFresnel")
fresnel.inputs["IOR"].default_value = n_water
nt.links.new(fresnel.outputs["Fac"], mix.inputs["Fac"])
nt.links.new(black.outputs["BSDF"], mix.inputs[1])    # fac=0 -> absorbed
nt.links.new(mirror.outputs["BSDF"], mix.inputs[2])   # fac=Fresnel -> reflect
nt.links.new(mix.outputs["Shader"], out.inputs["Surface"])
mesh_obj.data.materials.append(mat)

# --- uniform world background = I_sky -------------------------------------
world = bpy.data.worlds.new("flatsky")
scene.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
sky_rgb = tuple(spec.sky_rgb) if hasattr(spec, "sky_rgb") else (1.0, 1.0, 1.0)
bg.inputs["Color"].default_value = sky_rgb + (1.0,)
bg.inputs["Strength"].default_value = float(spec.I_sky)

# --- camera matching seapol._camera_rays ----------------------------------
# seapol's full horizontal FOV is 2*hfov_deg (rays span +/-tan(hfov_deg), i.e.
# +/-hfov_deg).  Set the focal length explicitly from the sensor width so the
# FOV is deterministic (cam_data.angle alone proved unreliable here).
cam_data = bpy.data.cameras.new("cam")
cam_data.sensor_fit = "HORIZONTAL"
cam_data.sensor_width = 36.0
cam_data.lens = (cam_data.sensor_width / 2.0) / float(np.tan(np.deg2rad(spec.hfov_deg)))
print(f"blender camera: lens={cam_data.lens:.1f}mm "
      f"-> hFOV={np.rad2deg(cam_data.angle_x):.2f} deg "
      f"(target {2*spec.hfov_deg:.2f} deg)")
cam_obj = bpy.data.objects.new("cam", cam_data)
scene.collection.objects.link(cam_obj)
# Aim the camera at the patch center with -Z along the view and Y up, the
# standard reliable Blender idiom (matches seapol's look/up basis).
origin = mathutils.Vector(cb["origin"])
direction = mathutils.Vector(cb["center"]) - origin
rot = direction.to_track_quat("-Z", "Y").to_matrix().to_4x4()
cam_obj.matrix_world = mathutils.Matrix.Translation(origin) @ rot
scene.camera = cam_obj

# --- render to linear EXR -------------------------------------------------
scene.render.resolution_x, scene.render.resolution_y = cb["W"], cb["H"]
scene.render.image_settings.file_format = "OPEN_EXR"
scene.render.image_settings.color_mode = "RGB"
scene.view_settings.view_transform = "Raw"   # no tonemap; keep linear radiance
scene.cycles.samples = 256
scene.render.filepath = out_exr
bpy.ops.render.render(write_still=True)

# Re-load the EXR through Blender's own image API and save a top-down RGB .npy
# next to it, so the comparison step needs no external EXR reader.  Blender's
# pixel buffer is row-major bottom-to-top RGBA; flip to top-down and drop alpha.
img = bpy.data.images.load(out_exr)
w, h = img.size
px = np.array(img.pixels[:], dtype=np.float32).reshape(h, w, img.channels)
px = px[::-1, :, :3]
np.save(out_exr[:-4] + ".npy" if out_exr.endswith(".exr") else out_exr + ".npy",
        px)
print(f"blender wrote {out_exr} and {out_exr[:-4]}.npy")
