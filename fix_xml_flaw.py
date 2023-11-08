'''
The XML files have been saved with a flaw
one of the attributes in saves as a type long 
it should be type double so the floating point 
values dont thow an error
'''
import xml.etree.ElementTree as ET
import os, glob

for ec in range(1,7):
    path = f'dataset/ec_{ec}/graph/'
    for filename in glob.glob(os.path.join(path,'*.gexf')):
        print(filename)
        #data = ET.parse(filename)
        #root = data.getroot()
        try:
            data = ET.parse(filename)
            root = data.getroot()
            #root[1][0][0].attrib['type'] = 'string'
            root[1][0][1].attrib['type'] = 'double'
            data.write(filename)
        except:
            print('error with:', filename)
