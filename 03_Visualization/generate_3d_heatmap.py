import slicer
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as colors
import json
import vtk
import math

# exec(open(r"C:/Users/matti/Documents/Vault/Exjobb/Code/Omics_Imaging_Integration/03_Visualization/generate_3d_heatmap.py").read())
all_features=["adrenal_gland_left",
"adrenal_gland_right",
"cardiac_fat",
"duodenum",
"esophagus",
"gallbladder",
"heart",
"inner_fat",
"intestine",
"kidney_left",
"kidney_right",
"liver",
"pancreas",
"sacrum",
"spleen",
"stomach",
"subcutaneous_fat",
"thyroid_gland",
"trachea",
"urinary_bladder",
"lung_upper_lobe_left",
"lung_lower_lobe_left",
"lung_upper_lobe_right",
"lung_middle_lobe_right",
"lung_lower_lobe_right",
"prostate",
"aorta",
"pulmonary_vein",
"brachiocephalic_trunk",
"subclavian_artery_right",
"subclavian_artery_left",
"common_carotid_artery_right",
"common_carotid_artery_left",
"brachiocephalic_vein_left",
"brachiocephalic_vein_right",
"atrial_appendage_left",
"superior_vena_cava",
"inferior_vena_cava",
"portal_vein_and_splenic_vein",
"iliac_artery_left",
"iliac_artery_right",
"iliac_vena_left",
"iliac_vena_right",
"humerus_left",
"humerus_right",
"scapula_left",
"scapula_right",
"clavicula_left",
"clavicula_right",
"femur_left",
"femur_right",
"hip_left",
"hip_right",
"spinal_cord",
"gluteus_maximus_left",
"gluteus_maximus_right",
"gluteus_medius_left",
"gluteus_medius_right",
"gluteus_minimus_left",
"gluteus_minimus_right",
"autochthon_left",
"autochthon_right",
"iliopsoas_left",
"iliopsoas_right",
"sternum",
"costal_cartilages",
"muscle",
"IVD",
"vertebra_body",
"vertebra_posterior_elements",
"spinal_channel",
"bone_other"]
csv_path = r"C:/Users/matti/Documents/Vault/Exjobb/results for heatmap.csv"
colormap = "Reds"

def get_model_from_folder(folder_name, model_name):
    shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
    folder_item = shNode.GetItemByName(folder_name)

    if folder_item == 0:
        print(f"Folder not found: {folder_name}")
        return None

    child_ids = vtk.vtkIdList()
    shNode.GetItemChildren(folder_item, child_ids, True)  # True = recursive

    for i in range(child_ids.GetNumberOfIds()):
        item_id = child_ids.GetId(i)
        node = shNode.GetItemDataNode(item_id)

        if node and node.IsA("vtkMRMLModelNode") and node.GetName() == model_name:
            return node

    print(f"Model not found in folder: {folder_name}/{model_name}")
    return None

def generate_heatmap(r2_col_male, r2_col_female, segmentation_name_male, segmentation_name_female, legend=False):
    seg_node_male = slicer.util.getNode(segmentation_name_male)
    seg_node_female = slicer.util.getNode(segmentation_name_female)

    df = pd.read_csv(csv_path)

    # Normalize color range
    vmin = 0
    vmax = 0.5
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(colormap)
    label_dict_path = "C:/Users/matti/Documents/Vault/Exjobb/Code/Omics_Imaging_Integration/03_Visualization/VIBE_labels.json"
    with open(label_dict_path) as f:
        label_dict = json.load(f)

    for seg_node, r2_col in zip([seg_node_male, seg_node_female], [r2_col_male, r2_col_female]):
        segmentation = seg_node.GetSegmentation()
        display_node = seg_node.GetDisplayNode()
        for _, row in df.iterrows():
            segment_name = row["feature"]
            if seg_node == seg_node_female and segment_name == "prostate":
                pass
            else:
                label = row["feature"]
                if segment_name in label_dict["name_to_id"]:
                    label_id = label_dict["name_to_id"][segment_name]
                    segment_name = "Segment_" + str(label_id)
                    segmentation.GetSegment(segment_name).SetName(label)
                value = row[r2_col]

                if pd.isna(value):
                    continue

                segment_id = segmentation.GetSegmentIdBySegmentName(label)



                if not segment_id:
                    print(f"Segment not found: {label}")
                    continue

                rgba = cmap(norm(value))
                rgb = rgba[:3]

                segmentation.GetSegment(segment_id).SetColor(rgb)
                display_node.SetSegmentOpacity(segment_id, 1.0)  # visible/opaque


    print("Done coloring segments by R2.")

    greenWidget = slicer.app.layoutManager().sliceWidget("Red")
    greenView = greenWidget.sliceView()
    renderer = greenView.renderWindow().GetRenderers().GetFirstRenderer()

    actors2D = renderer.GetActors2D()
    actors2D.InitTraversal()

    to_remove = []

    for i in range(actors2D.GetNumberOfItems()):
        actor = actors2D.GetNextActor2D()
        if actor and actor.IsA("vtkScalarBarActor"):
            to_remove.append(actor)

        if actor.IsA("vtkTextActor"):
            to_remove.append(actor)
            continue


    for actor in to_remove:
        renderer.RemoveActor2D(actor)



    greenView.renderWindow().Render()

    actors2D = renderer.GetActors2D()
    actors2D.InitTraversal()

    to_remove = []

    for i in range(actors2D.GetNumberOfItems()):
        actor = actors2D.GetNextActor2D()

        if actor and actor.IsA("vtkActor2D"):
            mapper = actor.GetMapper()
            if mapper:
                input_data = mapper.GetInput()
                if input_data and input_data.IsA("vtkPolyData"):
                    # crude filter: line-based actors (your ticks)
                    if input_data.GetNumberOfLines() > 0:
                        to_remove.append(actor)

    for actor in to_remove:
        renderer.RemoveActor2D(actor)

    greenView.renderWindow().Render()

    print("Removed legends.")

    if legend == True:

        lut = vtk.vtkLookupTable()
        lut.SetRange(vmin, vmax)
        lut.SetNumberOfTableValues(256)
        lut.Build()

        for i in range(256):
            x = i / 255
            x = x*vmax
            r, g, b, a = cmap(norm(x))
            lut.SetTableValue(i, r, g, b, a)

        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(lut)

        # labels every 0.1
        labels = vtk.vtkDoubleArray()
        tick_values = []


        start = math.ceil(vmin / 0.1) * 0.1
        end   = math.floor(vmax / 0.1) * 0.1

        x = start
        while x <= end + 1e-9:
            labels.InsertNextValue(round(x, 1))
            tick_values.append(round(x, 1))
            x += 0.1

        # scalar_bar.SetCustomLabels(labels)
        # scalar_bar.UseCustomLabelsOn()
        # scalar_bar.SetLabelFormat("%.1f")
        scalar_bar.DrawTickLabelsOff()
        labelText = scalar_bar.GetLabelTextProperty()
        scalar_bar.UnconstrainedFontSizeOn()

        greenWidget = slicer.app.layoutManager().sliceWidget("Red")
        greenView = greenWidget.sliceView()
        renderer = greenView.renderWindow().GetRenderers().GetFirstRenderer()

        renderer.AddActor2D(scalar_bar)
        scalar_bar.SetWidth(0.07)
        scalar_bar.SetHeight(0.55)
        scalar_bar.SetPosition(0.82, 0.23)
        scalar_bar.GetLabelTextProperty().SetFontSize(20)
        scalar_bar.SetLabelFormat("%.1f")
        greenView.renderWindow().Render()

        print("Added color legend.")
        greenView.renderWindow().Render()

        colorbar_tick_actors = []

        width, height = greenView.renderWindow().GetSize()

        bar_x = 0.82
        bar_y = 0.23
        bar_w = 0.07
        bar_h = 0.55

        # Empirical padding inside vtkScalarBarActor
        pad_bottom = 0.008
        pad_top = 0.001

        tick_len = 10  # pixels

        for val in tick_values:
            frac = (val - vmin) / (vmax - vmin)

            y_norm = bar_y + pad_bottom + frac * (bar_h - pad_bottom - pad_top)
            y = int(y_norm * height)

            x1 = int((bar_x + 0.027) * width)
            x2 = x1 + tick_len

            line = vtk.vtkLineSource()
            line.SetPoint1(x1, y, 0)
            line.SetPoint2(x2, y, 0)

            mapper = vtk.vtkPolyDataMapper2D()
            mapper.SetInputConnection(line.GetOutputPort())

            tick_actor = vtk.vtkActor2D()
            tick_actor.SetMapper(mapper)
            tick_actor.GetProperty().SetColor(0, 0, 0)
            tick_actor.GetProperty().SetLineWidth(2)

            renderer.AddActor2D(tick_actor)
            colorbar_tick_actors.append(tick_actor)

            text_actor = vtk.vtkTextActor()
            text_actor.SetInput(f"{val:.1f}")
            text_actor.GetTextProperty().SetFontSize(20)
            text_actor.GetTextProperty().SetColor(0, 0, 0)
            text_actor.GetTextProperty().ItalicOn()
            text_actor.GetTextProperty().BoldOn()

            text_actor.SetDisplayPosition(x2 + 5, y - 10)

            renderer.AddActor2D(text_actor)
            colorbar_tick_actors.append(text_actor)

        greenView.renderWindow().Render()


        print("Added color legend ticks.")
    feature = r2_col_male.replace("_m","")
    sliceViewName = "Red"
    filename = "C:/Users/matti/Downloads/" + feature + "_heatmap_" + colormap +".png"

    # Set view background to white
    view = slicer.app.layoutManager().sliceWidget(sliceViewName).sliceView()
    view.setBackgroundColor(qt.QColor.fromRgbF(1,1,1))
    view.forceRender()
    slicer.app.layoutManager().sliceWidget("Red").sliceView().setBackgroundColor(qt.QColor.fromRgbF(1,1,1))
    slicer.app.layoutManager().sliceWidget("Yellow").sliceView().setBackgroundColor(qt.QColor.fromRgbF(1,1,1))

    # Capture a screenshot
    import ScreenCapture
    cap = ScreenCapture.ScreenCaptureLogic()
    cap.captureImageFromView(view, filename)

def generate_heatmap_IMATs(r2_col_male, r2_col_female, segmentation_name_male, segmentation_name_female, legend=False):
    seg_node_male = slicer.util.getNode(segmentation_name_male)
    seg_node_female = slicer.util.getNode(segmentation_name_female)

    imat_features = ["autochthon_left_fat",
                     "autochthon_right_fat",
                     "gluteus_maximus_left_fat",
                     "gluteus_maximus_right_fat",
                     "gluteus_medius_left_fat",
                     "gluteus_medius_right_fat",
                     "gluteus_minimus_left_fat",
                     "gluteus_minimus_right_fat",
                     "iliopsoas_left_fat",
                     "iliopsoas_right_fat",
                     "muscle_fat",
                     "cardiac_fat"
                     ]

    df = pd.read_csv(csv_path)

    df = df[df["feature"].isin(imat_features)]

    mapping={"autochthon_left_fat": "autochthon_left",
           "autochthon_right_fat": "autochthon_right",
           "gluteus_maximus_left_fat": "gluteus_maximus_left",
           "gluteus_maximus_right_fat": "gluteus_maximus_right",
           "gluteus_medius_left_fat": "gluteus_medius_left",
           "gluteus_medius_right_fat": "gluteus_medius_right",
           "gluteus_minimus_left_fat": "gluteus_minimus_left",
           "gluteus_minimus_right_fat": "gluteus_minimus_right",
           "iliopsoas_left_fat": "iliopsoas_left",
           "iliopsoas_right_fat": "iliopsoas_right",
           "muscle_fat": "muscle",
             "cardiac_fat": "cardiac_fat"}

    df["feature"] = df["feature"].map(mapping)
    # Normalize color range
    vmin = 0
    vmax = 0.5
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(colormap)
    label_dict_path = "C:/Users/matti/Documents/Vault/Exjobb/Code/Omics_Imaging_Integration/03_Visualization/VIBE_labels.json"
    with open(label_dict_path) as f:
        label_dict = json.load(f)

    for seg_node, r2_col in zip([seg_node_male, seg_node_female], [r2_col_male, r2_col_female]):
        segmentation = seg_node.GetSegmentation()

        display_node = seg_node.GetDisplayNode()
        display_node.SetEdgeVisibility(True)
        display_node.SetLineWidth(5)

        for i in range(segmentation.GetNumberOfSegments()):
            segment_id = segmentation.GetNthSegmentID(i)
            display_node.SetSegmentOpacity(segment_id, 0)  # transparent

        for _, row in df.iterrows():
            segment_name = row["feature"]
            if seg_node == seg_node_female and segment_name == "prostate":
                pass
            else:
                label = row["feature"]
                if segment_name in label_dict["name_to_id"]:
                    label_id = label_dict["name_to_id"][segment_name]
                    segment_name = "Segment_" + str(label_id)
                    segmentation.GetSegment(segment_name).SetName(label)
                value = row[r2_col]

                if pd.isna(value):
                    continue

                segment_id = segmentation.GetSegmentIdBySegmentName(label)



                if not segment_id:
                    print(f"Segment not found: {label}")
                    continue

                rgba = cmap(norm(value))
                rgb = rgba[:3]

                segmentation.GetSegment(segment_id).SetColor(rgb)
                display_node.SetSegmentOpacity(segment_id, 1.0)  # visible/opaque
                display_node.SetSegmentOpacity3D(segment_id, 1.0)

    print("Done coloring segments by R2.")

    greenWidget = slicer.app.layoutManager().sliceWidget("Green")
    greenView = greenWidget.sliceView()
    renderer = greenView.renderWindow().GetRenderers().GetFirstRenderer()

    actors2D = renderer.GetActors2D()
    actors2D.InitTraversal()

    to_remove = []

    for i in range(actors2D.GetNumberOfItems()):
        actor = actors2D.GetNextActor2D()
        if actor and actor.IsA("vtkScalarBarActor"):
            to_remove.append(actor)

        if actor.IsA("vtkTextActor"):
            to_remove.append(actor)
            continue


    for actor in to_remove:
        renderer.RemoveActor2D(actor)



    greenView.renderWindow().Render()

    actors2D = renderer.GetActors2D()
    actors2D.InitTraversal()

    to_remove = []

    for i in range(actors2D.GetNumberOfItems()):
        actor = actors2D.GetNextActor2D()

        if actor and actor.IsA("vtkActor2D"):
            mapper = actor.GetMapper()
            if mapper:
                input_data = mapper.GetInput()
                if input_data and input_data.IsA("vtkPolyData"):
                    # crude filter: line-based actors (your ticks)
                    if input_data.GetNumberOfLines() > 0:
                        to_remove.append(actor)

    for actor in to_remove:
        renderer.RemoveActor2D(actor)

    greenView.renderWindow().Render()

    print("Removed legends.")

    if legend == True:

        lut = vtk.vtkLookupTable()
        lut.SetRange(vmin, vmax)
        lut.SetNumberOfTableValues(256)
        lut.Build()

        for i in range(256):
            x = i / 255
            x = x*vmax
            r, g, b, a = cmap(norm(x))
            lut.SetTableValue(i, r, g, b, a)

        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(lut)

        # labels every 0.1
        labels = vtk.vtkDoubleArray()
        tick_values = []


        start = math.ceil(vmin / 0.1) * 0.1
        end   = math.floor(vmax / 0.1) * 0.1

        x = start
        while x <= end + 1e-9:
            labels.InsertNextValue(round(x, 1))
            tick_values.append(round(x, 1))
            x += 0.1

        # scalar_bar.SetCustomLabels(labels)
        # scalar_bar.UseCustomLabelsOn()
        # scalar_bar.SetLabelFormat("%.1f")
        scalar_bar.DrawTickLabelsOff()
        labelText = scalar_bar.GetLabelTextProperty()
        scalar_bar.UnconstrainedFontSizeOn()

        greenWidget = slicer.app.layoutManager().sliceWidget("Green")
        greenView = greenWidget.sliceView()
        renderer = greenView.renderWindow().GetRenderers().GetFirstRenderer()

        renderer.AddActor2D(scalar_bar)
        scalar_bar.SetWidth(0.07)
        scalar_bar.SetHeight(0.55)
        scalar_bar.SetPosition(0.82, 0.23)
        scalar_bar.GetLabelTextProperty().SetFontSize(20)
        scalar_bar.SetLabelFormat("%.1f")
        greenView.renderWindow().Render()

        print("Added color legend.")
        greenView.renderWindow().Render()

        colorbar_tick_actors = []

        width, height = greenView.renderWindow().GetSize()

        bar_x = 0.82
        bar_y = 0.23
        bar_w = 0.07
        bar_h = 0.55

        # Empirical padding inside vtkScalarBarActor
        pad_bottom = 0.008
        pad_top = 0.001

        tick_len = 10  # pixels

        for val in tick_values:
            frac = (val - vmin) / (vmax - vmin)

            y_norm = bar_y + pad_bottom + frac * (bar_h - pad_bottom - pad_top)
            y = int(y_norm * height)

            x1 = int((bar_x + 0.027) * width)
            x2 = x1 + tick_len

            line = vtk.vtkLineSource()
            line.SetPoint1(x1, y, 0)
            line.SetPoint2(x2, y, 0)

            mapper = vtk.vtkPolyDataMapper2D()
            mapper.SetInputConnection(line.GetOutputPort())

            tick_actor = vtk.vtkActor2D()
            tick_actor.SetMapper(mapper)
            tick_actor.GetProperty().SetColor(0, 0, 0)
            tick_actor.GetProperty().SetLineWidth(2)

            renderer.AddActor2D(tick_actor)
            colorbar_tick_actors.append(tick_actor)

            text_actor = vtk.vtkTextActor()
            text_actor.SetInput(f"{val:.1f}")
            text_actor.GetTextProperty().SetFontSize(20)
            text_actor.GetTextProperty().SetColor(0, 0, 0)
            text_actor.GetTextProperty().ItalicOn()
            text_actor.GetTextProperty().BoldOn()

            text_actor.SetDisplayPosition(x2 + 5, y - 10)

            renderer.AddActor2D(text_actor)
            colorbar_tick_actors.append(text_actor)

        greenView.renderWindow().Render()


        print("Added color legend ticks.")
    feature = r2_col_male.replace("_m","")
    sliceViewName = "Green"
    filename = "C:/Users/matti/Downloads/" + feature + "_heatmap_" + colormap +"_IMAT.png"

    # Set view background to white
    view = slicer.app.layoutManager().sliceWidget(sliceViewName).sliceView()
    view.setBackgroundColor(qt.QColor.fromRgbF(1,1,1))
    view.forceRender()

    # Capture a screenshot
    import ScreenCapture
    cap = ScreenCapture.ScreenCaptureLogic()
    cap.captureImageFromView(view, filename)

def generate_heatmap_3D(interested_features, r2_col_male, r2_col_female, model_name_male, model_name_female, legend=False):
    for model_node in slicer.util.getNodesByClass("vtkMRMLModelNode"):
        display_node = model_node.GetDisplayNode()
        if display_node:
            display_node.SetVisibility(False)


    df = pd.read_csv(csv_path)
    df = df[df["feature"].isin(interested_features)]


    # Normalize color range
    vmin = 0
    vmax = 0.5
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(colormap)
    label_dict_path = "C:/Users/matti/Documents/Vault/Exjobb/Code/Omics_Imaging_Integration/03_Visualization/VIBE_labels.json"
    with open(label_dict_path) as f:
        label_dict = json.load(f)

    for model_node_name, r2_col in zip([model_name_male, model_name_female], [r2_col_male, r2_col_female]):

        for _, row in df.iterrows():
            segment_name = row["feature"]
            if model_node_name == model_name_female and segment_name == "prostate":
                pass
            else:

                label = row["feature"]
                if model_node_name == model_name_male:
                    gender_label = "_male"
                if model_node_name == model_name_female:
                    gender_label = "_female"
                label += gender_label
                value = row[r2_col]
                model_node = get_model_from_folder(model_node_name, label)
                if model_node is None:
                    print("No model for "+label)
                    continue
                display_node = model_node.GetDisplayNode()
                if display_node is None:
                    model_node.CreateDefaultDisplayNodes()
                    display_node = model_node.GetDisplayNode()

                rgba = cmap(norm(value))
                rgb = rgba[:3]
                display_node.SetColor(rgb)
                display_node.SetOpacity(1.0)
                display_node.SetVisibility(True)

        if "subcutaneous_fat" not in interested_features:
            model_node = get_model_from_folder(model_node_name, "subcutaneous_fat"+gender_label)
            display_node = model_node.GetDisplayNode()
            display_node.SetColor((0,0,0))
            display_node.SetOpacity(0.05)
            display_node.SetVisibility(True)

    print("Done coloring segments by R2.")

def generate_heatmap_3D_IMAT(r2_col_male, r2_col_female, model_name_male, model_name_female, legend=False):
    interested_features = ["autochthon_left_fat",
                     "autochthon_right_fat",
                     "gluteus_maximus_left_fat",
                     "gluteus_maximus_right_fat",
                     "gluteus_medius_left_fat",
                     "gluteus_medius_right_fat",
                     "gluteus_minimus_left_fat",
                     "gluteus_minimus_right_fat",
                     "iliopsoas_left_fat",
                     "iliopsoas_right_fat",
                     "cardiac_fat"
                     ]


    for model_node in slicer.util.getNodesByClass("vtkMRMLModelNode"):
        display_node = model_node.GetDisplayNode()
        if display_node:
            display_node.SetOpacity(0)


    df = pd.read_csv(csv_path)
    df = df[df["feature"].isin(interested_features)]
    mapping = {"autochthon_left_fat": "autochthon_left",
               "autochthon_right_fat": "autochthon_right",
               "gluteus_maximus_left_fat": "gluteus_maximus_left",
               "gluteus_maximus_right_fat": "gluteus_maximus_right",
               "gluteus_medius_left_fat": "gluteus_medius_left",
               "gluteus_medius_right_fat": "gluteus_medius_right",
               "gluteus_minimus_left_fat": "gluteus_minimus_left",
               "gluteus_minimus_right_fat": "gluteus_minimus_right",
               "iliopsoas_left_fat": "iliopsoas_left",
               "iliopsoas_right_fat": "iliopsoas_right",
               "cardiac_fat": "cardiac_fat"}

    df["feature"] = df["feature"].map(mapping)

    # Normalize color range
    vmin = 0
    vmax = 0.5
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(colormap)
    label_dict_path = "C:/Users/matti/Documents/Vault/Exjobb/Code/Omics_Imaging_Integration/03_Visualization/VIBE_labels.json"
    with open(label_dict_path) as f:
        label_dict = json.load(f)

    for model_node_name, r2_col in zip([model_name_male, model_name_female], [r2_col_male, r2_col_female]):

        for _, row in df.iterrows():
            segment_name = row["feature"]
            if model_node_name == model_name_female and segment_name == "prostate":
                pass
            else:
                label = row["feature"]
                value = row[r2_col]
                model_node = get_model_from_folder(model_node_name, label)
                if model_node is None:
                    print("No model for "+label)
                    continue
                display_node = model_node.GetDisplayNode()
                if display_node is None:
                    model_node.CreateDefaultDisplayNodes()
                    display_node = model_node.GetDisplayNode()

                rgba = cmap(norm(value))
                rgb = rgba[:3]
                display_node.SetColor(rgb)
                display_node.SetOpacity(1.0)
                display_node.SetVisibility(True)

        if "subcutaneous_fat" not in interested_features:
            model_node = get_model_from_folder(model_node_name, "subcutaneous_fat")
            display_node = model_node.GetDisplayNode()
            display_node.SetColor((0,0,0))
            display_node.SetOpacity(0.05)
            display_node.SetVisibility(True)

    print("Done coloring segments by R2.")
for measure, legend in zip(["ff_height", "vol_height", "ff", "vol"],[True, False, False, False]):
    generate_heatmap("prot_m_"+measure, "prot_f_"+measure, "2369783_processed_mask.nii.gz","4717924_processed_mask.nii.gz",legend=legend)
    generate_heatmap("met_m_"+measure, "met_f_"+measure, "2369783_processed_mask.nii.gz","4717924_processed_mask.nii.gz", legend=legend)


# for measure, legend in zip(["ff_height", "vol_height", "ff", "vol"],[True, False, False, False]):
#     generate_heatmap_IMATs("prot_m_"+measure, "prot_f_"+measure, "2369783_processed_mask.nii.gz","4717924_processed_mask.nii.gz",legend=legend)
#     generate_heatmap_IMATs("met_m_"+measure, "met_f_"+measure, "2369783_processed_mask.nii.gz","4717924_processed_mask.nii.gz", legend=legend)

bones = ["sacrum",
"humerus_left",
"humerus_right",
"scapula_left",
"scapula_right",
"clavicula_left",
"clavicula_right",
"femur_left",
"femur_right",
"hip_left",
"hip_right",
"spinal_cord",
"sternum",
"costal_cartilages",
"IVD",
"vertebra_body",
"vertebra_posterior_elements",
"spinal_channel",
"bone_other"
]

organs = ["adrenal_gland_left",
"adrenal_gland_right",
"gallbladder",
"kidney_left",
"kidney_right",
"liver",
"pancreas",
"spleen",
"thyroid_gland",
"urinary_bladder",
"prostate"
]

lungs=["trachea",
"lung_upper_lobe_left",
"lung_lower_lobe_left",
"lung_upper_lobe_right",
"lung_middle_lobe_right",
"lung_lower_lobe_right"
]

cardiovascular=["cardiac_fat",
"heart",
"aorta",
"pulmonary_vein",
"brachiocephalic_trunk",
"subclavian_artery_right",
"subclavian_artery_left",
"common_carotid_artery_right",
"common_carotid_artery_left",
"brachiocephalic_vein_left",
"brachiocephalic_vein_right",
"atrial_appendage_left",
"superior_vena_cava",
"inferior_vena_cava",
"portal_vein_and_splenic_vein",
"iliac_artery_left",
"iliac_artery_right",
"iliac_vena_left",
"iliac_vena_right"
]

fats = ["inner_fat",
"subcutaneous_fat"
]

digestion = ["duodenum",
"esophagus",
"intestine",
"stomach"
]

muscles = ["gluteus_maximus_left",
"gluteus_maximus_right",
"gluteus_medius_left",
"gluteus_medius_right",
"gluteus_minimus_left",
"gluteus_minimus_right",
"autochthon_left",
"autochthon_right",
"iliopsoas_left",
"iliopsoas_right",
"muscle"
]

muscle_fats = ["gluteus_maximus_left_fat",
"gluteus_maximus_right_fat",
"gluteus_medius_left_fat",
"gluteus_medius_right_fat",
"gluteus_minimus_left_fat",
"gluteus_minimus_right_fat",
"autochthon_left_fat",
"autochthon_right_fat",
"iliopsoas_left_fat",
"iliopsoas_right_fat",
"muscle_fat"]

# generate_heatmap_3D(cardiovascular,"prot_m_vol_height", "prot_f_vol_height", "2369783_processed_mask.nii.gz-models","4717924_processed_mask.nii.gz-models")
# generate_heatmap_3D_IMAT("prot_m_vol_height", "prot_f_vol_height", "2369783_processed_mask.nii.gz-models","4717924_processed_mask.nii.gz-models")


