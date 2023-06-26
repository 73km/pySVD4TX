Steps:

Get the pdb and *.map files ready and place them in "input_files" directory. The "input_files" directory needs to be created in the same directory where the svd.py program resides.

run find_coordinates.sh to find the orthogonal and fractional coordinates of desired residue or ligand for each map/pdb. This will generate a coordinate file with both orthogonal and fractional coordinates of each atom of the selected residue. There will also be a orthogonal coordinate file which will be used by the carving out program. Orthogonal coordinate file for all the subunits/residues of same map/pdb can be combined together

run mask_out.py to carve out the residue/ligand of interest. Check properly the input map file, the coordinate file and the subunit. Also set the radius of the mask around each atom. These are case sensitive and will be used by later steps.

If the unit cells remain unchanged in all the time delays, skip this step. Else, run carve_out_box.sh to carve out a box around the residue of interest. The program requires masked map from above step and positive fractional coordinates of an atom for reference.

run SVD.py to do SVD analysis. The input map array must be in chornologically ascending order. The program expects the maps in "input_files" directory. They should contain original maps, masked maps and/or carved maps.
