This Script creates a metabolic network based on an annotated reference network and a
provided network.

Input: A top-right-corner connection matrix of the network which should be compared to the reference:
Name: metabolic_network.csv

Row/column names or organs: 
Adrenal Glands
Bone Marrow	
Brain	
Colon	
Heart	
Kidney	
Liver	
Lung	
Muscle	
Pancreas	
Small Intestine	
Spleen	
Thyroid	
WAT

The provided network can also have only some of these organs.

Parameters:
threshold: Edge values above this value will be considered in the reference network building.
default = 0.3.

Output:
reference_network.html: The reference network.
metabolic_network_comparison.html: The comparison network between the provided and reference network.
metabolic_network_comparison.csv: The connection matrix of the comparison network.