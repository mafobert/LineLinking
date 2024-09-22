#! /usr/bin/env python

# This is my script to connect linear and curvilinear features in a feature class
# Coordinate system must be UTM
# shapefile must contain the following fields in the following order:
# FID, Shape, Cat, Azimuth, xstart, ystart, xend, yend

##### For CURVILINEAR - connect lines from start to end (add on, don't create new line) ######

#import arcpy
import math
arcpy.env.workspace = "C:/Users/mfobert/Desktop/Manic"

# assign values to link distance and azimuth thresholds (also patches)
d_link = 10000			#metres 
d_width = 3000			#(metres) orthogonal to the line
az_thres = 35			#degrees (+/-)
overlap = 12000			#overlap in metres between moving kernelx
k_width = 40000			#kernelx width
k_height = 40000			#kernelx height

# input dataset
input_fc= arcpy.env.workspace + "/AA_0001_final_7500_gt750_Spl_out_gt500.shp" #this is the dataset I am working on
print "input file: " + input_fc
fc_name="AA_0001_final_7500_gt750_Spl_out_gt500"
print "feature class name " + fc_name

# get list of field names in input file and save in fieldnames variable
# The Describe function returns a Describe Object, with multiple properties, such as datatype, fields, indexes, and many others.
dsc = arcpy.Describe(input_fc)
fields = dsc.fields
fieldnames = [field.name for field in fields]
print "These are the fieldnames of " + input_fc + ":  " 
print fieldnames
# Create outputx file based off of input file - has the same fieldList and spatial reference but no values
# CreateFeatureclass_management(out_path, out_name, {geometry_type}, {tempxlate;tempxlate...}, {has_m}, {has_z}, {spatial_reference}, {config_keyword}, {spatial_grid_1}, {spatial_grid_2}, {spatial_grid_3}
# creates feature class for lines
arcpy.CreateFeatureclass_management(arcpy.env.workspace, "outputx.shp", "POLYLINE", input_fc, "SAME_AS_TEMPLATE", "SAME_AS_TEMPLATE", input_fc)
# creates feature class for polygons so I can visually test the results
arcpy.CreateFeatureclass_management(arcpy.env.workspace, "corridorx.shp", "POLYGON", input_fc, "SAME_AS_TEMPLATE", "SAME_AS_TEMPLATE", input_fc)
# creates feature class for kernelx
arcpy.CreateFeatureclass_management(arcpy.env.workspace, "kernelx.shp", "POLYGON", "", "", "", input_fc)

result_fc = "outputx.shp"
poly_fc = "corridorx.shp"
kernelx_fc = "kernelx.shp"

# get input feature class' extent
left = dsc.extent.XMin
right = dsc.extent.XMax
bottom = dsc.extent.YMin
top = dsc.extent.YMax
print "These are the extents: "
print "left: %s \nright: %s \nbottom: %s \ntop: %s " % (left, right, bottom, top)

# column counter
g = 1
# row counter
h = 0
# overlap control in y direction
q = 0
ky3 = top
tempx_num = 0

# put kernal in a loop to move over the input data. 1:Upper Left, 2:Upper Right, 3:Lower Right, 4:Lower Left. Therefore: ky2 = ky1, kx3 = kx2, kx4 = kx1, ky4 = ky3
while  ky3 > (bottom - k_height/2):
   ky1 = ky3 + (q * overlap)
   ky3 = ky1 - k_height
   h += 1
   # overlap control in x direction
   p = 0
   q = 1
   kx2 = left
   while kx2 <= (right + k_width/2):
      print "Working on row:  " + str(h) + ",  patch:  " + str(g)
      if arcpy.Exists(arcpy.env.workspace + "/tempx_list" + str(tempx_num) + ".shp"):
         print "The tempx_list" + str(tempx_num) + "file exists and a new one will be created due to lock issue"
         tempx_num += 1

      # setup kernelx boundaries
      kx1 = kx2 - (p * overlap)
      kx2 = kx1 + k_width
      g += 1
      p = 1
      # build kernelx to pass over input feature extents
      kernelx_cords = [[kx1, ky1], [kx2, ky1], [kx2, ky3], [kx1, ky3], [kx1, ky1]]
      kernelx_array = arcpy.Array()	#create empty array object
      kernelx_point = arcpy.Point()	#create empty point object
      for k_coords in kernelx_cords:	#populate point object with X and Y coordinates from kernelx_cords and append to array object
         kernelx_point.X=k_coords[0]
         kernelx_point.Y=k_coords[1]
         kernelx_array.add(kernelx_point)
      kernelx_point=arcpy.Polygon(kernelx_array)	#creates polygon object using array object
      with arcpy.da.InsertCursor(kernelx_fc, ['SHAPE@']) as kernelx_Cursor:
         kernelx_Cursor.insertRow([kernelx_point])

      # Only select the lines that intersect the kernelx (polygon) - the results are not saved (they are only a temporary layer)
      try:
            batch = arcpy.SelectLayerByLocation_management(fc_name,'INTERSECT',"kernelx") #i think this may be the issue; kernelx is whole area, kernelx_fc is moving window
      except: 
            print "The error is at the batch step"
            print(arcpy.GetMessages())
      # put selected lines into a layer tempx_list (along with all their attributes) 
      # CopyFeatures_management(in_features, out_feature_class, {config_keyword}, {spatial_grid_1}, {spatial_grid_2}, {spatial_grid_3})
      arcpy.CopyFeatures_management(batch, 'tempx_list'+str(tempx_num))

      # create update cursor to save results to outputx file - use 'with' to keep from having the dataset locked
      print "You have made it to A"
      row_index = 0
      poly_Cursor = arcpy.da.InsertCursor(poly_fc, ['SHAPE@'])
      poly_CursorA = arcpy.da.InsertCursor(poly_fc, ['SHAPE@'])
      with arcpy.da.InsertCursor(result_fc, fieldnames) as inCursor:
         # read in parameters of the first row of the input dataset
         with arcpy.da.SearchCursor('tempx_list'+str(tempx_num), fieldnames) as sCursor:
            with arcpy.da.SearchCursor('tempx_list'+str(tempx_num), fieldnames) as sCursor2:
               for row in sCursor:
                  print "you have made it to B" #, this is the current row_index: " + str(row_index)
                  row_index += 1
                  FID1 = row[0]
                  cat = row[2]
                  az_1 = row[3]
                  xstart_1 = row[5]
                  ystart_1 = row[6]
                  xend_1 = row[7]
                  yend_1 = row[8]
                  # convert azimuth to 0 -> 180 deg (difference between GRASS GIS and ArcGIS)
                  if (az_1) >= 180: 
                     az_test = az_1 - 180
                  else:
                     az_test = az_1
                  # get azimuth range
                  az_max = az_test + az_thres
                  az_min_eval = az_test - az_thres
                  # Need to take into account boundary condition when az_min is negative
                  if (az_min_eval) < 0: 
                     az_min = 360 - (-1*az_min_eval)
                  else:
                     az_min = az_min_eval
                  # find northern node and assign as X2,Y2 - southern node will be X1,Y1
                  if (ystart_1 < yend_1):
                     X1 = xstart_1
                     Y1 = ystart_1
                     X2 = xend_1
                     Y2 = yend_1
                  else:
                     X1 = xend_1
                     Y1 = yend_1
                     X2 = xstart_1
                     Y2 = ystart_1
                  # convert degrees to radians to use in cos and sin functions
                  az_rad = (math.pi/180)*az_test
                  # set-up corridorx (create polygon around line)
                  if (az_test <= 90):
                     x1mid = X1 - (d_link * math.sin(az_rad))
                     y1mid = Y1 - (d_link * math.cos(az_rad))
                     x2mid = X2 + (d_link * math.sin(az_rad))
                     y2mid = Y2 + (d_link * math.cos(az_rad))
                     px1 = x1mid + ((d_width/2)*math.sin((math.pi/2)-az_rad))
                     py1 = y1mid - ((d_width/2)*math.cos((math.pi/2)-az_rad))
                     px2 = x2mid + ((d_width/2)*math.sin((math.pi/2)-az_rad))
                     py2 = y2mid - ((d_width/2)*math.cos((math.pi/2)-az_rad))
                     px3 = x2mid - ((d_width/2)*math.sin((math.pi/2)-az_rad))
                     py3 = y2mid + ((d_width/2)*math.cos((math.pi/2)-az_rad))
                     px4 = x1mid - ((d_width/2)*math.sin((math.pi/2)-az_rad))
                     py4 = y1mid + ((d_width/2)*math.cos((math.pi/2)-az_rad))
                     list_cords = [ [px1, py1], [px2, py2], [px3, py3], [px4, py4], [px1, py1]]
                  elif (az_test > 90 and az_test <= 180):
                     x1mid = X1 + d_link*(math.sin((math.pi)-az_rad))
                     y1mid = Y1 - d_link*(math.cos((math.pi)-az_rad))
                     x2mid = X2 - d_link*(math.sin((math.pi)-az_rad))
                     y2mid = Y2 + d_link*(math.cos((math.pi)-az_rad))
                     px1 = x1mid + (d_width/2)*(math.sin(az_rad-(math.pi/2)))
                     py1 = y1mid + (d_width/2)*(math.cos(az_rad-(math.pi/2)))
                     px2 = x2mid + (d_width/2)*(math.sin(az_rad-(math.pi/2)))
                     py2 = y2mid + (d_width/2)*(math.cos(az_rad-(math.pi/2)))
                     px3 = x2mid - (d_width/2)*(math.sin(az_rad-(math.pi/2)))
                     py3 = y2mid - (d_width/2)*(math.cos(az_rad-(math.pi/2)))
                     px4 = x1mid - (d_width/2)*(math.sin(az_rad-(math.pi/2)))
                     py4 = y1mid - (d_width/2)*(math.cos(az_rad-(math.pi/2)))
                     list_cords = [ [px1, py1], [px2, py2], [px3, py3], [px4, py4], [px1, py1]]
                  array = arcpy.Array()		#create empty array object
                  point_obj = arcpy.Point()	#create empty point object
                  for coords in list_cords:	#populate point object with X and Y coordinates from list_cords and append to array object
                     point_obj.X=coords[0]
                     point_obj.Y=coords[1]
                     array.add(point_obj)
                  poly_obj=arcpy.Polygon(array)	#creates polygon object using array object
                  poly_Cursor.insertRow([poly_obj])
                  # now set-up perpendicular box (polygon) around line
                  # Perpendicular is to test for stacking
                  if (az_test > 0 and az_test < 90):
                     beta = (math.pi/2) - az_rad
                     A1X1 = xstart_1 - (d_width*3)*(math.cos(beta))
                     A1X2 = xstart_1 + (d_width*3)*(math.cos(beta))
                     A1Y1 = ystart_1 - (d_width*3)*(math.sin(beta))
                     A1Y2 = ystart_1 + (d_width*3)*(math.sin(beta))
                     A2X1 = xend_1 - (d_width*3)*(math.cos(beta))
                     A2X2 = xend_1 + (d_width*3)*(math.cos(beta))
                     A2Y1 = yend_1 - (d_width*3)*(math.sin(beta))
                     A2Y2 = yend_1 + (d_width*3)*(math.sin(beta))
                     list_cordsA1 = [ [A1X1, A1Y1], [A2X1, A2Y1], [A2X2, A2Y2], [A1X2, A1Y2]]
                  elif (az_test == 0 or az_test == 180):
                     A1X1 = xstart_1 - (d_width*3)
                     A1X2 = xstart_1 + (d_width*3)
                     A1Y1 = ystart_1
                     A1Y2 = ystart_1
                     A2X1 = xend_1 - (d_width*3)
                     A2X2 = xend_1 + (d_width*3)
                     A2Y1 = yend_1
                     A2Y2 = yend_1
                     list_cordsA1 = [ [A1X1, A1Y1], [A2X1, A2Y1], [A2X2, A2Y2], [A1X2, A1Y2]]
                  elif (az_test == 90):
                     A1X1 = xstart_1
                     A1X2 = xstart_1
                     A1Y1 = ystart_1 - (d_width*3)
                     A1Y2 = ystart_1 + (d_width*3)
                     A2X1 = xend_1
                     A2X2 = xend_1
                     A2Y1 = yend_1 - (d_width*3)
                     A2Y2 = yend_1 + (d_width*3)
                     list_cordsA1 = [ [A1X1, A1Y1], [A2X1, A2Y1], [A2X2, A2Y2], [A1X2, A1Y2]]
                  elif (az_test > 90 and az_test < 180):
                     beta = (math.pi) - az_rad
                     A1X1 = xstart_1 - (d_width*3)*(math.cos(beta))
                     A1X2 = xstart_1 + (d_width*3)*(math.cos(beta))
                     A1Y1 = ystart_1 - (d_width*3)*(math.sin(beta))
                     A1Y2 = ystart_1 + (d_width*3)*(math.sin(beta))
                     A2X1 = xend_1 - (d_width*3)*(math.cos(beta))
                     A2X2 = xend_1 + (d_width*3)*(math.cos(beta))
                     A2Y1 = yend_1 - (d_width*3)*(math.sin(beta))
                     A2Y2 = yend_1 + (d_width*3)*(math.sin(beta))
                     list_cordsA1 = [ [A1X1, A1Y1], [A2X1, A2Y1], [A2X2, A2Y2], [A1X2, A1Y2]]
                  arrayA1 = arcpy.Array()	#create empty array object
                  point_objA1 = arcpy.Point()	#create empty point object
                  for coords in list_cordsA1:	#populate point object with X and Y coordinates from list_cords and append to array object
                     point_objA1.X=coords[0]
                     point_objA1.Y=coords[1]
                     arrayA1.add(point_objA1)
                  poly_objA1=arcpy.Polygon(arrayA1)	#creates polygon object using array object
                  poly_CursorA.insertRow([poly_objA1])
                  # At this point, the first line has been read in and a polygon has been created around it, and well
                  # as an additional one to check for stacking
                  # Now, read in parameters of FID+1 for comparison
                  # advance sCursor2 to the current row in sCursor
                  sCursor2.reset()        
                  for i in xrange(0, row_index):
                  sCursor2.next() 
                  # loop over all next rows
                  for row2 in sCursor2:
                     FID2 = row2[0]
                     az_n = row2[3]
                     xstart_n = row2[5]
                     ystart_n = row2[6]
                     xend_n = row2[7]
                     yend_n = row2[8]
                     # convert azimuth to 0 -> 180 deg (difference between GRASS GIS and ArcGIS)
                     if (az_n) >= 180: 
                        az_n2 = az_n - 180
                     else:
                        az_n2 = az_n
                     # another boundary condition for when az_1 is close to 0 (i.e., 358)
                     az_nBound = az_n2 + 180
                     # criteria for connecting features
                     # first check intersection of nodes with polygon, then azimuth threshold
                     if FID1 != FID2:	#this is to ensure you don't compare the same row to itself
                        # Create points list array
                        pt_C = arcpy.Point(xstart_n, ystart_n)
                        pt_D = arcpy.Point(xend_n, yend_n)
                        # Now check if one or both of the nodes are within the polygon. if yes, continue 
                        if (poly_obj.contains(pt_C)) or (poly_obj.contains(pt_D)):	
                        # Either node is within polygon, NOW check for azimuth threshold
                           if ((az_n2 <= az_max) and (az_n2 >= az_min)) or ((az_min_eval < 0) and (az_nBound <= 360) and (az_nBound >= az_min)) or ((az_min_eval < 0) and (az_n2 <= az_max) and (az_n2 >= 0)) or ((az_nBound <= az_max) and (az_nBound >= az_min)):
                              #Azimuth threshold met, continue
                              az = 2
                              # Here is where I will check for stacking (does the perpendicular polygon/rectangle intersect with either of the points)? then I will check for 'stepping'
                              if (poly_objA1.contains(pt_C)) or (poly_objA1.contains(pt_D)):
                              # Line _n stacks (slightly parallel) with line _1, so skip
                                 pass		#pass if lines are at all parallel
                              else:
                                 # Check which point is both within the polygon AND is closest to _1
                                 # get the minimum distances
                                 check_AC = math.sqrt(math.pow((xstart_n - xstart_1),2) + math.pow((ystart_n - ystart_1),2))
                                 check_AD = math.sqrt(math.pow((xend_n - xstart_1),2) + math.pow((yend_n - ystart_1),2))
                                 check_BC = math.sqrt(math.pow((xstart_n - xend_1),2) + math.pow((ystart_n - yend_1),2))
                                 check_BD = math.sqrt(math.pow((xend_n - xend_1),2) + math.pow((yend_n - yend_1),2))
                                 check_min = min (check_AC, check_AD, check_BC, check_BD)
                                 if (poly_obj.contains(pt_C)) and (min(check_AC, check_BC) < min(check_AD, check_BD)): #check if this point is closest to _1
                                    # point C is in the polygon AND is the closest to _1, proceed, else check point D
                                    # Now I need to know if the point closest to C is A or B
                                    if (check_AC == check_min):
                                       # This let's us know point C is closest to A 
                                       # if check_min == 0, the two nodes are coexisting
                                       if (check_AC == 0):
                                          # just add back in the two lines, give a value of 10
                                          print "Nodes A and C are coincident, just keep both lines"
                                          inCursor.insertRow((1, 1, cat, az, 10, xstart_1, ystart_1, xend_1, yend_1))
                                          inCursor.insertRow((1, 1, cat, az, 10, xstart_n, ystart_n, xend_n, yend_n))
                                       else:
                                          # Now we need to see if the new (connected) line would be 'stepped'
                                          # this is done by calculating the angle between vector AA and AC
                                          unEditedAz = (math.pi/180)*az_1
                                          hypot = 10 
                                          AextendX = xstart_1 - hypot*(math.sin(unEditedAz))
                                          AextendY = ystart_1 - hypot*(math.cos(unEditedAz))
                                          # now we create vectors AAextend and AC (I extended AA but did not need to)
                                          AAvectX = AextendX - xstart_1
                                          AAvectY = AextendY - ystart_1
                                          ACvectX = xstart_n - xstart_1
                                          ACvectY = ystart_n - ystart_1
                                          # find dot product of AA dot AC
                                          dot_prod = (AAvectX*ACvectX) + (AAvectY*ACvectY)
                                          AA_mag = math.sqrt((AAvectX*AAvectX)+(AAvectY*AAvectY))
                                          AC_mag = math.sqrt((ACvectX*ACvectX)+(ACvectY*ACvectY))
                                          value_test = dot_prod/(AA_mag*AC_mag)
                                          if (value_test < -1) or (value_test > 1):
                                             print "ERROR HERE at the AA_mag part - outside arccos domain! "
                                             pass
                                          # get magnitudes
                                          else:
                                             # equation is rho = arccos(AA dot AC/magnitude(AA)*magnitude(AC))
                                             rho = math.acos(value_test)
                                             # convert rho to degrees
                                             rho_check = (180*rho)/math.pi
                                             # check if rho is within the threshold for non-stepping
                                             if (rho_check <= az_thres/6):
                                                #yes, continue
                                                xstart_new = xstart_1
                                                ystart_new = ystart_1
                                                xend_new = xstart_n
                                                yend_new = ystart_n
                                                print "you have added in connection AC"
                                                inCursor.insertRow((1, 1, cat, az, 1, xstart_1, ystart_1, xend_1, yend_1))
                                                inCursor.insertRow((1, 1, cat, az, 1, xstart_n, ystart_n, xend_n, yend_n))
                                                inCursor.insertRow((1, 1, cat, az, 6, xstart_new, ystart_new, xend_new, yend_new))
                                             else:
                                                print "Skipping as lines are stepped"
                                                pass	# new line would be stepped, so skip
                                    elif (check_BC == check_min):
                                       # this let's us know point C is closest to B 
                                       # if check_min == 0, the two nodes are coexisting
                                       if (check_BC == 0):
                                          # just keep the two lines
                                          print "Nodes B and C are coincident, just keep both lines"
                                          inCursor.insertRow((1, 1, cat, az, 11, xstart_1, ystart_1, xend_1, yend_1))
                                          inCursor.insertRow((1, 1, cat, az, 11, xstart_n, ystart_n, xend_n, yend_n))
                                       else:
                                          # Now we need to see if the new (connected) line would be 'stepped'
                                          # this is done by calculating the angle between vector BB and BC
                                          unEditedAz = (math.pi/180)*az_1
                                          hypot = 10 
                                          BextendX = xend_1 + hypot*(math.sin(unEditedAz))
                                          BextendY = yend_1 + hypot*(math.cos(unEditedAz))
                                          # now we create vectors BBextend and BC (I extended BB but did not need to)
                                          BBvectX = BextendX - xend_1
                                          BBvectY = BextendY - yend_1
                                          BCvectX = xstart_n - xend_1
                                          BCvectY = ystart_n - yend_1
                                          # find dot product of AA dot AC
                                          dot_prod = (BBvectX*BCvectX) + (BBvectY*BCvectY)
                                          # get magnitudes
                                          BB_mag = math.sqrt((BBvectX*BBvectX)+(BBvectY*BBvectY))
                                          BC_mag = math.sqrt((BCvectX*BCvectX)+(BCvectY*BCvectY))
                                          value_test = dot_prod/(BB_mag*BC_mag)
                                          if (value_test < -1) or (value_test > 1):
                                             print "ERROR at the BC_mag step ! outside arccos domain!!!!!!!!!!"
                                             pass
                                          else:
                                             # equation is rho = arccos(BB dot BC/magnitude(BB)*magnitude(BC))
                                             rho = math.acos(value_test)
                                             # convert rho to degrees
                                             rho_check = (180*rho)/math.pi
                                             # check if rho is within the threshold for non-stepping
                                             if (rho_check <= az_thres/6):
                                                #yes, continue
                                                xstart_new = xend_1
                                                ystart_new = yend_1
                                                xend_new = xstart_n
                                                yend_new = ystart_n
                                                print "you have added in connection BC"
                                                inCursor.insertRow((1, 1, cat, az, 2, xstart_1, ystart_1, xend_1, yend_1))
                                                inCursor.insertRow((1, 1, cat, az, 2, xstart_n, ystart_n, xend_n, yend_n))
                                                inCursor.insertRow((1, 1, cat, az, 6, xstart_new, ystart_new, xend_new, yend_new))
                                             else:
                                                print "Skipping as lines are stepped"
                                                pass		# new line would be stepped, so skip
                                    else:
                                       print "point C is in the polygon AND is the closest to _1, but there is an error"
                                 elif (poly_obj.contains(pt_D)) and (min(check_AD, check_BD) <= min(check_AC, check_BC)):
                                    # point D is in the polygon AND is the closest to _1, proceed, else post error
                                    # Now I need to know if the point closest to D is A or B
                                    if (check_AD == check_min):
                                       # This let's us know point D is closest to A 
                                       # if check_min == 0, the two nodes are coexisting
                                       if (check_AD == 0):
                                          print "Nodes A and D are coincident, keep both lines"
                                          inCursor.insertRow((1, 1, cat, az, 9, xstart_1, ystart_1, xend_1, yend_1))
                                          inCursor.insertRow((1, 1, cat, az, 9, xstart_n, ystart_n, xend_n, yend_n))
                                       else:
                                          # this let's us know point D is closest to A 
                                          # Now we need to see if the new (connected) line would be 'stepped'
                                          # this is done by calculating the angle between vector AA and AD
                                          unEditedAz = (math.pi/180)*az_1
                                          hypot = 10 
                                          AextendX = xstart_1 - hypot*(math.sin(unEditedAz))
                                          AextendY = ystart_1 - hypot*(math.cos(unEditedAz))
                                          # now we create vectors AAextend and AD (I extended AA but did not need to)
                                          AAvectX = AextendX - xstart_1
                                          AAvectY = AextendY - ystart_1
                                          ADvectX = xend_n - xstart_1
                                          ADvectY = yend_n - ystart_1
                                          # find dot product of AA dot AD
                                          dot_prod = (AAvectX*ADvectX) + (AAvectY*ADvectY)
                                          # get magnitudes
                                          AA_mag = math.sqrt((AAvectX*AAvectX)+(AAvectY*AAvectY))
                                          AD_mag = math.sqrt((ADvectX*ADvectX)+(ADvectY*ADvectY))
                                          value_test = dot_prod/(AA_mag*AD_mag)
                                          if (value_test < -1) or (value_test > 1):
                                             print " ERROR HERE AD_mag step !! outside arccos domain!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                                             pass
                                          else:
                                             # equation is rho = arccos(AA dot AC/magnitude(AA)*magnitude(AC))
                                             rho = math.acos(value_test)
                                             # convert rho to degrees
                                             rho_check = (180*rho)/math.pi
                                             # check if rho is within the threshold for non-stepping
                                             if (rho_check <= az_thres/6):
                                                #yes, continue
                                                xstart_new = xstart_1
                                                ystart_new = ystart_1
                                                xend_new = xend_n
                                                yend_new = yend_n
                                                print "you have added in connection AD"
                                                inCursor.insertRow((1, 1, cat, az, 3, xstart_1, ystart_1, xend_1, yend_1))
                                                inCursor.insertRow((1, 1, cat, az, 3, xstart_n, ystart_n, xend_n, yend_n))
                                                inCursor.insertRow((1, 1, cat, az, 6, xstart_new, ystart_new, xend_new, yend_new))
                                             else:
                                                print "Skipping as lines are stepped"
                                                pass		# new line would be stepped, so skip
                                    elif (check_BD == check_min):
                                       # this let's us know point D is closest to B 
                                       # if check_min == 0, the two nodes are coexisting
                                       if (check_BD == 0):
                                          print "Nodes B and D are coincident, keep both lines"
                                          inCursor.insertRow((1, 1, cat, az, 10, xstart_1, ystart_1, xend_1, yend_1))
                                          inCursor.insertRow((1, 1, cat, az, 10, xstart_n, ystart_n, xend_n, yend_n))
                                       else:
                                          # Now we need to see if the new (connected) line would be 'stepped'
                                          # this is done by calculating the angle between vector BB and BD
                                          unEditedAz = (math.pi/180)*az_1
                                          hypot = 10 
                                          BextendX = xend_1 + hypot*(math.sin(unEditedAz))
                                          BextendY = yend_1 + hypot*(math.cos(unEditedAz))
                                          # now we create vectors BBextend and BC (I extended BB but did not need to)
                                          BBvectX = BextendX - xend_1
                                          BBvectY = BextendY - yend_1
                                          BDvectX = xend_n - xend_1
                                          BDvectY = yend_n - yend_1
                                          # find dot product of AA dot AC
                                          dot_prod = (BBvectX*BDvectX) + (BBvectY*BDvectY)
                                          # get magnitudes
                                          BB_mag = math.sqrt((BBvectX*BBvectX)+(BBvectY*BBvectY))
                                          BD_mag = math.sqrt((BDvectX*BDvectX)+(BDvectY*BDvectY))
                                          value_test = dot_prod/(BB_mag*BD_mag)
                                          if (value_test < -1) or (value_test > 1):
                                             print "ERROR at BD_mag  !!! outside arccos domain!!"
                                             pass
                                          else:
                                             # equation is rho = arccos(BB dot BD/magnitude(BB)*magnitude(BD))
                                             rho = math.acos(value_test)
                                             # convert rho to degrees
                                             rho_check = (180*rho)/math.pi
                                             # check if rho is within the threshold for non-stepping
                                             if (rho_check <= az_thres/6):
                                                #yes, continue
                                                xstart_new = xend_1
                                                ystart_new = yend_1
                                                xend_new = xend_n
                                                yend_new = yend_n
                                                print "you have added in connection BD"
                                                inCursor.insertRow((1, 1, cat, az, 4, xstart_1, ystart_1, xend_1, yend_1))
                                                inCursor.insertRow((1, 1, cat, az, 4, xstart_n, ystart_n, xend_n, yend_n))
                                                inCursor.insertRow((1, 1, cat, az, 6, xstart_new, ystart_new, xend_new, yend_new))
                                             else:
                                                print "Skipping as lines are stepped"
                                                pass		# new line would be stepped, so skip
                                    else:
                                       print "point D is in the polygon AND is the closest to _1, but there is an error"
                                 else:
                                    print "either node detected as being inside, but niether observed to be closest: Thus, closest node not within the polygon"
                           else:		#pass if azimuth not within threshold range
                              pass
                                                  else:			#pass if neither of the nodes are within the polygon
                          pass			
							else:				#pass if FID1 = FID2
								pass
		arcpy.Delete_management("kernelx.shp")
		del kernelx_Cursor
		arcpy.CreateFeatureclass_management(arcpy.env.workspace, "kernelx.shp", "POLYGON", "", "", "", input_fc)

print "Now removing duplicates"
#copied and edited for my dataset from https://community.esri.com/thread/163776
table_rows = []  
with arcpy.da.UpdateCursor("outputx", fieldnames) as dcursor:
	for drow in dcursor:  
		#print "this is the xstart value:   " + str(drow[4])
		if drow[4:] in table_rows:  
			print "Deleting record: FID = {}".format(drow[0])  
			dcursor.deleteRow()  
		else:  
			table_rows.append(drow[4:])  
del table_rows

print "Converting outputx to shapefile"
#ArcToolbox > Data Management Tools > Features > XY To Line
# Set local variables
input_table = arcpy.env.workspace + "/outputx.shp"
out_lines = arcpy.env.workspace + "/AA_00001_Oct.shp"
#XY To Line
arcpy.XYToLine_management(input_table,out_lines,"xstart","ystart","xend","yend","GEODESIC","cat") #try fid instead of cat
print "Done"

