
inmuscle = {}

inmuscle["2200_Wrists_and_hand_bones_cortical"]=[]
inmuscle["3700_Ankles_and_foot_cortical"]=[]
inmuscle["4500_Scapulae_cortical"]=[]
inmuscle["10000_Lymphatic_nodes_ET"]=[]
inmuscle["10100_Lymphatic_nodes_thoracic"]=[]
inmuscle["10200_Lymphatic_nodes_head"]=[]
inmuscle["10300_Lymphatic_nodes_trunk"]=[]
inmuscle["10400_Lymphatic_nodes_arms"]=[]
inmuscle["10500_Lymphatic_nodes_legs"]=[]


inRSTlist=[         

"100_Adrenal_left",
"200_Adrenal_right",
"303_ET1_surface",
"405_ET2_surface",
"500_Oral_mucosa_tongue",
###"501_Oral_mucosa_mouth_floor",
###"600_Oral_mucosa_lip_and_cheeks",
"700_Trachea",
"808_BB1_surface",
"900_Blood_in_large_arteries",
"910_Blood_in_large_veins",
"1300_Humeri_cortical",#, "1500_Humeri_medullary_cavity"
"1900_Ulnae_and_radii_cortical",
"2200_Wrists_and_hand_bones_cortical",
"2400_Clavicles_cortical",
"2600_Cranium_cortical",
###"2600_Cranium_cortical_in"
###"2600_Cranium_cortical_surrounding_frontal_sinus"
"2800_Femora_cortical",
"3400_Tibiae_fibulae_and_patellae_cortical",
###"3700_Ankles_and_foot_cortical",
"3900_Mandible_cortical",
"4100_Pelvis_cortical",
"4300_Ribs_cortical",
	#"4500_Scapulae_cortical",
"4700_Cervical_spine_cortical",
"4900_Thoracic_spine_cortical",
"5100_Lumbar_spine_cortical",
"5300_Sacrum_cortical",
"5500_Sternum_cortical",
"5700_Cartilage_costal",
"5800_Cartilage_discs",
"6100_Brain",
"6200_Breast_left_adipose_tissue",
"6300_Breast_left_glandular_tissue",
"6400_Breast_right_adipose_tissue",
"6500_Breast_right_glandular_tissue",
"7000_Gall_bladder_wall",
"7203_Stomach_wall_surface",
"7403_Small_intestine_wall_surface",


"8700_Heart_wall",
###"8800_Blood_in_heart_chamber",
"8900_Kidney_left_cortex",
"9200_Kidney_right_cortex",
"9500_Liver",

"9700_Lung_(AI)_left",
"9900_Lung_(AI)_right",

"10000_Lymphatic_nodes_ET",
"10100_Lymphatic_nodes_thoracic",
"10200_Lymphatic_nodes_head",
"10300_Lymphatic_nodes_trunk",
"10400_Lymphatic_nodes_arms",
"10500_Lymphatic_nodes_legs",
"10600_Muscle",
"11002_Oesophagus_wall_surface",
"11100_Ovary_left",
"11200_Ovary_right",
"11300_Pancreas",
"11400_Pituitary_gland",
"12000_Salivary_glands_left",
"12100_Salivary_glands_right",
"12600_Spinal_cord",
"12700_Spleen",
###"12801_Teeth_retention_region",
"12800_Teeth",
"13100_Thymus",
"13200_Thyroid",
"13300_Tongue_upper_(food)",
"13301_Tongue_lower_surface",
"13400_Tonsils",
"13500_Ureter_left",
"13600_Ureter_right",
"13700_Urinary_bladder_wall",
"13900_Uterus",

"3700_Ankles_and_foot_cortical",

"7602_Ascending_colon_wall_surface",
"7802_Transverse_colon_wall_right_surface",
"8002_Transverse_colon_wall_left_surface",
"8202_Descending_colon_wall_surface",
"8402_Sigmoid_colon_wall_surface",
"8600_Rectum_wall",

	
	
]

#inRSTlist=[
#	"2200_Wrists_and_hand_bones_cortical",
#	"10600_Muscle"
#]

inRST={}
for rst in inRSTlist:
    inRST[rst]=[]



list_objs_names = {

    "100_Adrenal_left": {"out":"100_Adrenal_left","in": {}},
    "200_Adrenal_right": {"out":"200_Adrenal_right","in": {}},

    "300_ET1_8um": {"out":"300_ET1_8um", "in": {"14000_ET1_contents_0um_(air)":[]}},
    "301_ET1_40um": {"out":"301_ET1_40um","in": {"300_ET1_8um":[]}},
    "302_ET1_50um": {"out":"302_ET1_50um","in": {"301_ET1_40um":[]}},
    "303_ET1_surface": {"out":"303_ET1_surface","in": {"302_ET1_50um":[]}},

    "400_ET2_0um": {"out":"400_ET2_0um","in":{"14000_ET2_contents_-15um_(air)":[]}},
    "401_ET2_40um": {"out":"401_ET2_40um","in":{"400_ET2_0um":[]}},
    "402_ET2_50um": {"out":"402_ET2_50um","in":{"401_ET2_40um":[]}},
    "403_ET2_55um": {"out":"403_ET2_55um","in":{"402_ET2_50um":[]}},
    "404_ET2_65um": {"out":"404_ET2_65um","in":{"403_ET2_55um":[]}},
    "405_ET2_surface": {"out":"405_ET2_surface","in":{"404_ET2_65um":[]}},

    "500_Oral_mucosa_tongue": {"out":"500_Oral_mucosa_tongue","in":{"13301_Tongue_lower_-200um":[]}},
    "501_Oral_mucosa_mounth_floor": {"out":"501_Oral_mucosa_mouth_floor","in":{}},
    "600_Oral_mucosa_lip_and_cheeks": {"out":"600_Oral_mucosa_lip_and_cheeks","in":{}},
    "700_Trachea": {"out":"700_Trachea","in":{"14000_Trachea_contents_(air)":[]}},

    "800_BB1_-6um": {"out":"800_BB1_-6um","in":{"14000_BB1_contents_-11um_(air)":[]}},
    "801_BB1_0um": {"out":"801_BB1_0um","in":{"800_BB1_-6um":[]}},
    "802_BB1_10um": {"out":"802_BB1_10um","in":{"801_BB1_0um":[]}},
    "803_BB1_35um": {"out":"803_BB1_35um","in":{"802_BB1_10um":[]}},
    "804_BB1_40um": {"out":"804_BB1_40um","in":{"803_BB1_35um":[]}},
    "805_BB1_50um": {"out":"805_BB1_50um","in":{"804_BB1_40um":[]}},
    "806_BB1_60um": {"out":"806_BB1_60um","in":{"805_BB1_50um":[]}},
    "807_BB1_70um": {"out":"807_BB1_70um","in":{"806_BB1_60um":[]}},
    "808_BB1_surface": {"out":"808_BB1_surface","in":{"807_BB1_70um":[]}},

    "900_Blood_in_large_arteries": {"out":"900_Blood_in_large_arteries", "in":{}},
    "910_Blood_in_large_veins": {"out":"910_Blood_in_large_veins","in":{}},

    "1300_Humeri_cortical": {"out":"1300_Humeri_cortical","in":{"1700_Humeri_lower_spongiosa":[],"1500_Humeri_medullary_cavity":[],"1400_Humeri_upper_spongiosa":[]}},
    "1400_Humeri_upper_spongiosa": {"out":"1400_Humeri_upper_spongiosa","in":{}},
    "1500_Humeri_medullary_cavity": {"out":"1500_Humeri_medullary_cavity","in":{}},
    "1700_Humeri_lower_spongiosa": {"out":"1700_Humeri_lower_spongiosa","in":{}},

    "1900_Ulnae_and_radii_cortical": {"out":"1900_Ulnae_and_radii_cortical","in":{"2000_Ulnae_and_radii_spongiosa":[],"2100_Ulnae_and_radii_medullary_cavity":[]}},
    "2000_Ulnae_and_radii_spongiosa": {"out":"2000_Ulnae_and_radii_spongiosa","in":{}},
    "2100_Ulnae_and_radii_medullary_cavity": {"out":"2100_Ulnae_and_radii_medullary_cavity","in":{}},
    "2200_Wrists_and_hand_bones_cortical": {"out":"2200_Wrists_and_hand_bones_cortical","in":{"2300_Wrists_and_hand_bones_spongiosa":[]}},
    "2300_Wrists_and_hand_bones_spongiosa": {"out":"2300_Wrists_and_hand_bones_spongiosa","in":{}},
    "2400_Clavicles_cortical": {"out":"2400_Clavicles_cortical","in":{"2500_Clavicles_spongiosa":[]}},
    "2500_Clavicles_spongiosa": {"out":"2500_Clavicles_spongiosa","in":{}},
    

    "2600_Cranium_cortical_in": {"out":"2600_Cranium_cortical_in","in":{"11600_RST_in":[]}},
    "2600_Cranium_cortical_surrounding_frontal_sinus": {"out":"2600_Cranium_cortical_surrounding_frontal_sinus","in":{"14000_Air_remaining":[]}},
    "2600_Cranium_cortical": { "out": "2600_Cranium_cortical","in": {"2700_Cranium_spongiosa":[]}},
    "2700_Cranium_spongiosa": {"out":"2700_Cranium_spongiosa","in":{"2600_Cranium_cortical_in":[],"2600_Cranium_cortical_surrounding_frontal_sinus":[]}},

    "2800_Femora_cortical": {"out":"2800_Femora_cortical","in":{"2900_Femora_upper_spongiosa":[],"3000_Femora_medullary_cavity":[],"3200_Femora_lower_spongiosa":[]}},
    "2900_Femora_upper_spongiosa": {"out":"2900_Femora_upper_spongiosa","in":{}},
    "3000_Femora_medullary_cavity": {"out":"3000_Femora_medullary_cavity","in":{}},  
    "3200_Femora_lower_spongiosa": {"out":"3200_Femora_lower_spongiosa","in":{}},
    "3400_Tibiae_fibulae_and_patellae_cortical": {"out":"3400_Tibiae_fibulae_and_patellae_cortical","in":{"3500_Tibiae_fibulae_and_patellae_spongiosa":[],"3600_Tibiae_fibulae_and_patellae_medullary_cavity":[]}},
    "3500_Tibiae_fibulae_and_patellae_spongiosa": {"out":"3500_Tibiae_fibulae_and_patellae_spongiosa","in":{}},
    "3600_Tibiae_fibulae_and_patellae_medullary_cavity": {"out":"3600_Tibiae_fibulae_and_patellae_medullary_cavity","in":{}},
    "3700_Ankles_and_foot_cortical": {"out":"3700_Ankles_and_foot_cortical","in":{"3800_Ankles_and_foot_spongiosa":[]}},
    "3800_Ankles_and_foot_spongiosa": {"out":"3800_Ankles_and_foot_spongiosa","in":{}},
    "3900_Mandible_cortical": {"out":"3900_Mandible_cortical","in":{"4000_Mandible_spongiosa":[]}},
    "4000_Mandible_spongiosa": {"out":"4000_Mandible_spongiosa","in":{}},
    "4100_Pelvis_cortical": {"out":"4100_Pelvis_cortical","in":{"4200_Pelvis_spongiosa":[]}},
    "4200_Pelvis_spongiosa": {"out":"4200_Pelvis_spongiosa","in":{}},
    "4300_Ribs_cortical": {"out":"4300_Ribs_cortical","in":{"4400_Ribs_spongiosa":[]}},
    "4400_Ribs_spongiosa": {"out":"4400_Ribs_spongiosa","in":{}},
    "4500_Scapulae_cortical": {"out":"4500_Scapulae_cortical","in":{"4600_Scapulae_spongiosa":[]}},
    "4600_Scapulae_spongiosa": {"out":"4600_Scapulae_spongiosa","in":{}},
    "4700_Cervical_spine_cortical": {"out":"4700_Cervical_spine_cortical","in":{"4800_Cervical_spine_spongiosa":[]}},
    "4800_Cervical_spine_spongiosa": {"out":"4800_Cervical_spine_spongiosa","in":{}},
    "4900_Thoracic_spine_cortical": {"out":"4900_Thoracic_spine_cortical","in":{"5000_Thoracic_spine_spongiosa":[]}},
    "5000_Thoracic_spine_spongiosa": {"out":"5000_Thoracic_spine_spongiosa","in":{}},
    "5100_Lumbar_spine_cortical": {"out":"5100_Lumbar_spine_cortical","in":{"5200_Lumbar_spine_spongiosa":[]}},
    "5200_Lumbar_spine_spongiosa": {"out":"5200_Lumbar_spine_spongiosa","in":{}},
    "5300_Sacrum_cortical": {"out":"5300_Sacrum_cortical","in":{"5400_Sacrum_spongiosa":[]}},
    "5400_Sacrum_spongiosa": {"out":"5400_Sacrum_spongiosa","in":{}},
    "5500_Sternum_cortical": {"out":"5500_Sternum_cortical","in":{"5600_Sternum_spongiosa":[]}},
    "5600_Sternum_spongiosa": {"out":"5600_Sternum_spongiosa","in":{}},
    "5700_Cartilage_costal": {"out":"5700_Cartilage_costal","in":{}},
    "5800_Cartilage_discs": {"out":"5800_Cartilage_discs","in":{}},

    "6100_Brain": {"out":"6100_Brain","in":{}},
    "6200_Breast_left_adipose_tissue": {"out":"6200_Breast_left_adipose_tissue","in":{}},
    "6300_Breast_left_glandular_tissue": {"out":"6300_Breast_left_glandular_tissue","in":{}},
    "6400_Breast_right_adipose_tissue": {"out":"6400_Breast_right_adipose_tissue","in":{}},
    "6500_Breast_right_glandular_tissue": {"out":"6500_Breast_right_glandular_tissue","in":{}},

    "6600_Eye_lens_sensitive_left": {"out":"6600_Eye_lens_sensitive_left","in":{}},
    "6601_Eye_lens_insensitive_left": {"out":"6601_Eye_lens_insensitive_left","in":{}},
    "6700_Cornea_left": {"out":"6700_Cornea_left","in":{"6701_Aqueous_left":[],"6600_Eye_lens_sensitive_left":[],"6601_Eye_lens_insensitive_left":[],"6702_Vitreous_left":[]}},
    "6701_Aqueous_left": {"out":"6701_Aqueous_left","in":{}},
    "6702_Vitreous_left": {"out":"6702_Vitreous_left","in":{}},
    "6800_Eye_lens_sensitive_right": {"out":"6800_Eye_lens_sensitive_right","in":{}},
    "6801_Eye_lens_insensitive_right": {"out":"6801_Eye_lens_insensitive_right","in":{}},
    "6900_Cornea_right": {"out":"6900_Cornea_right","in":{"6800_Eye_lens_sensitive_right":[],"6801_Eye_lens_insensitive_right":[],"6901_Aqueous_right":[],"6902_Vitreous_right":[]}},
    "6901_Aqueous_right": {"out":"6901_Aqueous_right","in":{}},
    "6902_Vitreous_right": {"out":"6902_Vitreous_right","in":{}},

    "7000_Gall_bladder_wall": {"out":"7000_Gall_bladder_wall","in":{"7100_Gall_bladder_contents":[]}},
    "7100_Gall_bladder_contents": {"out":"7100_Gall_bladder_contents","in":{}},

    "7200_Stomach_wall_60um": {"out":"7200_Stomach_wall_60um","in":{"7300_Stomach_contents":[]}},
    "7201_Stomach_wall_100um": {"out":"7201_Stomach_wall_100um","in":{"7200_Stomach_wall_60um":[]}},
    "7202_Stomach_wall_300um": {"out":"7202_Stomach_wall_300um","in":{"7201_Stomach_wall_100um":[]}},
    "7203_Stomach_wall_surface": {"out":"7203_Stomach_wall_surface","in":{"7202_Stomach_wall_300um":[]}},
    "7300_Stomach_contents": {"out":"7300_Stomach_contents","in":{}},

    "7400_Small_intestine_wall_130um": {"out":"7400_Small_intestine_wall_130um","in":{"7500_Small_intestine_contents_0um":[]}},
    "7401_Small_intestine_wall_150um": {"out":"7401_Small_intestine_wall_150um","in":{"7400_Small_intestine_wall_130um":[]}},
    "7402_Small_intestine_wall_200um": {"out":"7402_Small_intestine_wall_200um","in":{"7401_Small_intestine_wall_150um":[]}},
    "7403_Small_intestine_wall_surface": {"out":"7403_Small_intestine_wall_surface","in":{"7402_Small_intestine_wall_200um":[]}},
    "7500_Small_intestine_contents_0um": {"out":"7500_Small_intestine_contents_0um","in":{"7501_Small_intestine_contents_-500um":[]}},
    "7501_Small_intestine_contents_-500um": {"out":"7501_Small_intestine_contents_-500um","in":{}},

    "7600_Ascending_colon_wall_280um": {"out":"7600_Ascending_colon_wall_280um","in":{"7700_Ascending_colon_contents":[]}},
    "7601_Ascending_colon_wall_300um": {"out":"7601_Ascending_colon_wall_300um","in":{"7600_Ascending_colon_wall_280um":[]}},
    "7602_Ascending_colon_wall_surface": {"out":"7602_Ascending_colon_wall_surface","in":{"7601_Ascending_colon_wall_300um":[]}},
    "7700_Ascending_colon_contents": {"out":"7700_Ascending_colon_contents","in":{}},

    "7800_Transverse_colon_wall_right_280um": {"out":"7800_Transverse_colon_wall_right_280um","in":{"7900_Transverse_colon_contents_right":[]}},
    "7801_Transverse_colon_wall_right_300um": {"out":"7801_Transverse_colon_wall_right_300um","in":{"7800_Transverse_colon_wall_right_280um":[]}},
    "7802_Transverse_colon_wall_right_surface": {"out":"7802_Transverse_colon_wall_right_surface","in":{"7801_Transverse_colon_wall_right_300um":[]}},
    "7900_Transverse_colon_contents_right": {"out":"7900_Transverse_colon_contents_right","in":{}},
    "8000_Transverse_colon_wall_left_280um": {"out":"8000_Transverse_colon_wall_left_280um","in":{"8100_Transverse_colon_contents_left":[]}},
    "8001_Transverse_colon_wall_left_300um": {"out":"8001_Transverse_colon_wall_left_300um","in":{"8000_Transverse_colon_wall_left_280um":[]}},
    "8002_Transverse_colon_wall_left_surface": {"out":"8002_Transverse_colon_wall_left_surface","in":{"8001_Transverse_colon_wall_left_300um":[]}},
    "8100_Transverse_colon_contents_left": {"out":"8100_Transverse_colon_contents_left","in":{}},

    "8200_Descending_colon_wall_280um": {"out":"8200_Descending_colon_wall_280um","in":{"8300_Descending_colon_contents":[]}},
    "8201_Descending_colon_wall_300um": {"out":"8201_Descending_colon_wall_300um","in":{"8200_Descending_colon_wall_280um":[]}},
    "8202_Descending_colon_wall_surface": {"out":"8202_Descending_colon_wall_surface","in":{"8201_Descending_colon_wall_300um":[]}},
    "8300_Descending_colon_contents": {"out":"8300_Descending_colon_contents","in":{}},
    
    "8400_Sigmoid_colon_wall_280um": {"out":"8400_Sigmoid_colon_wall_280um","in":{"8500_Sigmoid_colon_contents":[]}},
    "8401_Sigmoid_colon_wall_300um": {"out":"8401_Sigmoid_colon_wall_300um","in":{"8400_Sigmoid_colon_wall_280um":[]}},
    "8402_Sigmoid_colon_wall_surface": {"out":"8402_Sigmoid_colon_wall_surface","in":{"8401_Sigmoid_colon_wall_300um":[]}},
    "8500_Sigmoid_colon_contents": {"out":"8500_Sigmoid_colon_contents","in":{}},

    "8600_Rectum_wall": {"out":"8600_Rectum_wall","in":{}},
    "8700_Heart_wall": {"out":"8700_Heart_wall","in":{"8800_Blood_in_heart_chamber":[]}},
    "8800_Blood_in_heart_chamber": {"out":"8800_Blood_in_heart_chamber","in":{}},
    "8900_Kidney_left_cortex": {"out":"8900_Kidney_left_cortex","in":{"9000_Kidney_left_medulla":[]}},
    "9000_Kidney_left_medulla": {"out":"9000_Kidney_left_medulla","in":{"9100_Kidney_left_pelvis":[]}},
    "9100_Kidney_left_pelvis": {"out":"9100_Kidney_left_pelvis","in":{}},
    "9200_Kidney_right_cortex": {"out":"9200_Kidney_right_cortex","in":{"9300_Kidney_right_medulla":[]}},
    "9300_Kidney_right_medulla": {"out":"9300_Kidney_right_medulla","in":{"9400_Kidney_right_pelvis":[]}},
    "9400_Kidney_right_pelvis": {"out":"9400_Kidney_right_pelvis","in":{}},
    "9500_Liver": {"out":"9500_Liver","in":{}},
    "9700_Lung_AI_left": {"out":"9700_Lung_(AI)_left","in":{}},
    "9900_Lung_AI_right": {"out":"9900_Lung_(AI)_right","in":{}},

    "10000_Lymphatic_nodes_ET": {"out":"10000_Lymphatic_nodes_ET","in":{}},
    "10100_Lymphatic_nodes_thoracic": {"out":"10100_Lymphatic_nodes_thoracic","in":{}},
    "10200_Lymphatic_nodes_head": {"out":"10200_Lymphatic_nodes_head","in":{}},
    "10300_Lymphatic_nodes_trunk": {"out":"10300_Lymphatic_nodes_trunk","in":{}},
    "10400_Lymphatic_nodes_arms": {"out":"10400_Lymphatic_nodes_arms","in":{}},
    "10500_Lymphatic_nodes_legs": {"out":"10500_Lymphatic_nodes_legs","in":{}},

    "10600_Muscle": {"out":"10600_Muscle","in":inmuscle},

    "11000_Oesophagus_wall_190um": {"out":"11000_Oesophagus_wall_190um","in":{"11003_Oesophagus_wall_contents":[]}},
    "11001_Oesophagus_wall_200um": {"out":"11001_Oesophagus_wall_200um","in":{"11000_Oesophagus_wall_190um":[]}},
    "11002_Oesophagus_wall_surface": {"out":"11002_Oesophagus_wall_surface","in":{"11001_Oesophagus_wall_200um":[]}},
    "11003_Oesophagus_wall_contents": {"out":"11003_Oesophagus_wall_contents","in":{}},

    "11100_Ovary_left": {"out":"11100_Ovary_left","in":{}},
    "11200_Ovary_right": {"out":"11200_Ovary_right","in":{}},
    "11300_Pancreas": {"out":"11300_Pancreas","in":{}},
    "11400_Pituitary_gland": {"out":"11400_Pituitary_gland","in":{}},     

    "11600_RST": {"out":"11600_RST","in":inRST},
    "11600_RST_in": {"out":"11600_RST_in","in":{}},              

    "12000_Salivary_glands_left": {"out":"12000_Salivary_glands_left","in":{}},
    "12100_Salivary_glands_right": {"out":"12100_Salivary_glands_right","in":{}},

    "12200_Skin_sensitive": {"out":"12200_Skin_100um","in":{"12201_Skin_50um":[]}},
    "12201_Skin_insensitive": {"out":"12201_Skin_50um","in":{"11600_RST":[]}},

    "12600_Spinal_cord": {"out":"12600_Spinal_cord","in":{}},
    "12700_Spleen": {"out":"12700_Spleen","in":{}},
    "12800_Teeth": {"out":"12800_Teeth","in":{}},
    "12801_Teeth_retention_region": {"out":"12801_Teeth_retention_region","in":{}},
    "13100_Thymus": {"out":"13100_Thymus","in":{}},
    "13200_Thyroid": {"out":"13200_Thyroid","in":{}},

    "13300_Tongue_upper_food": {"out":"13300_Tongue_upper_(food)","in":{}},
    "13301_Tongue_lower_surface": {"out":"13301_Tongue_lower_surface","in":{}},
    "13400_Tonsils": {"out":"13400_Tonsils","in":{}},

    "13500_Ureter_left": {"out":"13500_Ureter_left","in":{}},
    "13600_Ureter_right": {"out":"13600_Ureter_right","in":{}},

    "13700_Urinary_bladder_wall_insensitive": {"out":"13700_Urinary_bladder_wall","in":{"13701_Urinary_bladder_wall_185um":[]}},
    "13701_Urinary_bladder_wall_sensitive": {"out":"13701_Urinary_bladder_wall_185um","in":{"13700_Urinary_bladder_wall_116um":[]}},
    "13702_Urinary_bladder_wall_in": {"out":"13700_Urinary_bladder_wall_116um","in":{"13800_Urinary_bladder_contents":[]}},
    "13800_Urinary_bladder_contents": {"out":"13800_Urinary_bladder_contents","in":{}},
    "13900_Uterus": {"out":"13900_Uterus","in":{}},
    
    "14000_Air": {"out":"14000_Air_remaining","in":{}},
    "14000_BB1_contents_-11um_air": {"out":"14000_BB1_contents_-11um_(air)","in":{}},
    "14000_ET1_contents_0um_air": {"out":"14000_ET1_contents_0um_(air)","in":{}},
    "14000_ET2_contents_-15um_air": {"out":"14000_ET2_contents_-15um_(air)","in":{}},
    "14000_Trachea_contents_air": {"out":"14000_Trachea_contents_(air)","in":{}}


}

'''

    ,"12000_Salivary_glands_left": ["12000_Salivary_glands_left"]
    ,"12100_Salivary_glands_right": ["12100_Salivary_glands_right"]

    ,"12200_Skin_sensitive": ["12200_Skin_100um","12201_Skin_50um"]
    ,"12201_Skin_insensitive": ["12201_Skin_50um","11600_RST"]


    ,"12600_Spinal_cord": ["12600_Spinal_cord"]
    ,"12700_Spleen": ["12700_Spleen"]
    ,"12800_Teeth": ["12800_Teeth"]
    ,"12801_Teeth_retention_region": ["12801_Teeth_retention_region"]
    ,"13100_Thymus": ["13100_Thymus"]
    ,"13200_Thyroid": ["13200_Thyroid"]
    ,"13300_Tongue_upper_food": ["13300_Tongue_upper_(food)"]
    #,"13301_Tongue_lower_-200um": ["13301_Tongue_lower_-200um"]
    ,"13301_Tongue_lower_surface": ["13301_Tongue_lower_surface"]
    ,"13400_Tonsils": ["13400_Tonsils"]
    ,"13500_Ureter_left": ["13500_Ureter_left"]
    ,"13600_Ureter_right": ["13600_Ureter_right"]
    ,"13700_Urinary_bladder_wall_insensitive": ["13700_Urinary_bladder_wall","13701_Urinary_bladder_wall_185um"]
    ,"13701_Urinary_bladder_wall_sensitive": ["13701_Urinary_bladder_wall_185um","13700_Urinary_bladder_wall_116um"]
    ,"13702_Urinary_bladder_wall_in": ["13700_Urinary_bladder_wall_116um","13800_Urinary_bladder_contents"]
    ,"13800_Urinary_bladder_contents": ["13800_Urinary_bladder_contents"]
    ,"13900_Uterus": ["13900_Uterus"]
    ,"14000_Air": ["14000_Air_remaining"]
    ,"14000_BB1_contents_-11um_air": ["14000_BB1_contents_-11um_(air)"]
    ,"14000_ET1_contents_0um_air": ["14000_ET1_contents_0um_(air)"]
    ,"14000_ET2_contents_-15um_air": ["14000_ET2_contents_-15um_(air)"]
    ,"14000_Trachea_contents_air": ["14000_Trachea_contents_(air)"]
    #,"12202_Skin_insensitive_out": ["12200_Skin_100um","12200_Skin_surface"]
    }
'''