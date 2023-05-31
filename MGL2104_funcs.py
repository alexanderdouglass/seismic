# Bunch of useful functions for analyzing MGL2104 Data, specifically for work in analyzing the impacts of different seabed layers.

import numpy as np
from matplotlib import pyplot as plt
import scipy.ndimage as spnd
import binascii
import struct

ci = complex(0,1)

def readSegd_2104(filename):
    '''
    
    '''
    # Initializing variables
    class header:
        file_number = 0
        format_code = 0
        gen_consts = 0
        year = 0
        additional_header_blocks = 1
        julian_day = 0
        hour = 0
        minute = 0
        second = 0
        manufacturer_code = 0
        manufacturer_serial_number = 0
        bytes_per_scan = 0
        base_scan_interval = 0
        polarity_code = 0
        scans_in_block_exponent = 0
        scans_in_block_base = 0 
        record_type_code = 0
        record_length = 0
        scans_type_per_record = 0
        channel_set_per_scan_type = 0
        added_skew_fields = 0
        extended_header_length = 0
        external_header_length = 0
        SEGD_revision = 0
        general_trailers = 0
        channel_Set1_channels = 0
        channel_Set2_channels = 0

    # This version only supports reading files with format code 8058. Also this version only reads General Header #1 and General Header #2
    # It does not read any information from any other general headers. It skips reading skew headers, external header blocks, and extended
    # header

    # Mention the location of the file
    with open(filename, 'rb') as f:
        content = f.read()
        temp = binascii.hexlify(content)

        ####### Part 1 - Reading Header Data 

        # Initialising and Computing Variables
        header.file_number = (temp[0:4]).decode('utf-8')
        header.format_code = int(temp[4:8])
        if header.format_code != 8058:
            raise Exception('This File can only read SegD files having format Code 8058')
        header.gen_consts = temp[8:20]
        header.year = int('20' + temp[20:22].decode('utf-8'))
        header.additional_header_blocks = int(temp[22:23])
        header.julian_day = int(temp[23:26])
        header.hour = int(temp[26:28])
        header.minute = int(temp[28:30])
        header.second = int(temp[30:32])
        header.manufacturer_code = int(temp[32:34])
        header.manufacturer_serial_number = int(temp[34:38])
        header.bytes_per_scan = int(temp[38:44])
        header.base_scan_interval = int(temp[44:46])
        header.polarity_code = bin(int(temp[46:47].decode('utf-8')))
        header.record_type_code = bin(int(temp[50:51].decode('utf-8')))
        header.record_length = temp[51:54].decode('utf-8')
        header.scans_type_per_record = int(temp[54:56])
        header.channel_set_per_scan_type = temp[56:58].decode('utf-8')
        header.added_skew_fields = int(temp[58:60])
        header.extended_header_length = temp[60:62].decode('utf-8')
        header.external_header_length = temp[62:64].decode('utf-8') 

        if header.scans_type_per_record > 1:
            raise Exception('This version of readSegD only handles Seg-D files with a single scan type per record')

        # Reading General Header Block #2 
        if header.external_header_length == 'ff':
            header.external_header_length = int(temp[78:82].decode('utf-8'), 16)

        # if header.record_length == 'fff':
        #     header.record_length = int(temp[92:98])

        header.SEGD_revision = float(temp[84:88]) 
        header.general_trailers = temp[88:92].decode('utf-8')


        # This version does not read additional header blocks

        # Reading Scan Type Header
        count = 128 + (header.additional_header_blocks) * 32

        header.channel_Set1_channels = int(temp[count+16:count+20].decode('utf-8'))

        count = count + 64 * int(header.channel_set_per_scan_type)


        # Read Skew Fields
        if header.added_skew_fields > 0:
            count = count + header.added_skew_fields * 32 * 2
            print('Skipping Skew Fields')

        # Read Extended Header
        count = count + int(header.extended_header_length) * 32 * 2
        # print('Skipping Extended Header Blocks')


        # Read External Length
        count = count + int(header.external_header_length) * 32 *2
        # print('Skipping External Header Blocks')

        ###### Part 2 - Reading Trace Data 
        # count = count + 18
        # search = int(temp[count:count + 2].decode('utf-8'))
        
        count = count + 54
        samples_per_trace = (int(temp[count: count + 6], 16))
        # print(samples_per_trace)
        count = count + 434
        # print(count)
        
        add1 = 64
        add2 = 0*8*6001
#         print(temp[count + add1*0 + add2*0:count + add1*1 + add2*0])
#         print(temp[count + add1*1 + add2*1:count + add1*2 + add2*1])
#         print(temp[count + add1*2 + add2*2:count + add1*3 + add2*2])
#         print(temp[count + add1*3 + add2*3:count + add1*4 + add2*3])
        
        data1 = np.zeros((samples_per_trace, header.channel_Set1_channels))
        for j in range(0, header.channel_Set1_channels):
            for i in range(0, samples_per_trace):
                temp1 = temp[count:count + 8].decode('utf-8')
                data1[i][j] = struct.unpack('!f', bytes.fromhex(temp1))[0]
                count = count + 8
                # if j == 1:
                    # print(temp1)
                    # print(data1[i][j])
            count = count + 1*40 + 7*64

        return header, data1
            
        

def readP190_2104(filename):
    '''
    Read P190 Navigation file for MGL2104
    '''
    import re

    # filename = '/media/asd21/My Passport/MGL2104_PD12_CSs/MGL2104PD12.0.p190'

    RECEIVER_NUMBER = 1200

    class nav:
        numCables = 0
        shotStart = 0
        shotEnd = 0
        depth = []
        vesselX = []
        vesselY = []
        sourceX = []
        sourceY = []
        tailX = []
        tailY = []       
        receiverX = []
        receiverY = []
        receiverZ = []

    file = open(filename, 'r')
    lines = file.readlines()
    receiverX = []
    receiverY = []
    receiverZ = []
    for line in lines:
        if line.startswith('H0101GENERAL'):
            nav.numCables = int(re.sub(r'.*(\d)\sCABLE.*', r'\1 ', line))

        elif line.startswith('H2600Line'):
            nav.shotStart, nav.shotEnd = [int(i) for i in line.split() if i.isdigit()]

        elif line.startswith('VGL') or line.startswith('VMGL'):
            res = line.split('W')
            xyd = (re.sub(r'\.(\d)', r'.\1 ', res[1])).split()
            nav.vesselX.append(float(xyd[0]))
            nav.vesselY.append(float(xyd[1]))
            if len(xyd) >= 3:
                nav.depth.append(float(xyd[2]))
            else:
                nav.depth.append(0)

        elif line.startswith('SGL') or line.startswith('SMGL'):
            res = line.split('W')
            xyd = (re.sub(r'\.(\d)', r'.\1 ', res[1])).split()
            nav.sourceX.append(float(xyd[0]))
            nav.sourceY.append(float(xyd[1]))

        elif line.startswith('CGL') or line.startswith('CMGL'):
            res = line.split('W')
            xyd = (re.sub(r'\.(\d)', r'.\1 ', res[1])).split()
            nav.tailX.append(float(xyd[0]))
            nav.tailY.append(float(xyd[1]))

        elif line.startswith('R'):
          # Add the current receiver arrays to nav and clear them 
            if len(receiverX) == RECEIVER_NUMBER * 1: #nav.numCables:
                nav.receiverX.append(receiverX)
                nav.receiverY.append(receiverY)
                nav.receiverZ.append(receiverZ)
                receiverX = []
                receiverY = []
                receiverZ = []

            length = 26
            receivers = [line[i: i + length] for i in range(1, length * 3 + 1, length)]
            for rec in receivers:
                xyz = (re.sub(r'\.(\d)', r'.\1 ', rec[4:])).split()
                receiverX.append(float(xyz[0]))
                receiverY.append(float(xyz[1]))
                if len(xyz) == 3:
                    receiverZ.append(float(xyz[2]))
                else:
                    receiverZ.append(0)

    return nav




def SigExp(geometry, volumes, angle, ghost=False):
    '''
    This function generates a source signal from an airgun array based on individual airgun signal estimates and propagation angle.
    INPUTS
        -geometry: [x,y,z], list of 3 numpy arrays describing airgun geometry.
        -volumes: numpy array of airgun array sizes
        -angle: propagation direction from airgun array, in degrees
        -rtmp: 
        -ghost (optional): Whether or not to include the surface reflection 
    OUTPUTS
        -srcSig: Source signal composed of airgun signal linear combination
    '''

    rtmp = 1e4
    AirgunSigs = np.load('AirgunSigsMGL2104_New.npz')
    AirgunDict = {}
    for volume in [40,60,90,120,180,220,360]:
        AirgunDict['Airgun'+str(volume)] = AirgunSigs['Airgun'+str(volume)]

    c = 1500
    Fs = 500
    Tend = 1
    df = 1/Tend
    
    f = np.arange(0,Fs+df,df)
    tmpLen = len(f)
    if np.mod(tmpLen,2) == 0:
        f = np.append(np.arange(0,Fs/2,df),np.arange(-Fs/2,0,df))
    else:
        f = np.append(np.arange(0,Fs/2+df,df),np.arange(-Fs/2,0,df))

    p = np.zeros(len(f),dtype=complex)

    xg = geometry[0]
    yg = geometry[1]
    zg = geometry[2]

    tshift = 2*rtmp/c-0.15

    for ii in range(len(volumes)):

        tmp = np.copy(AirgunDict['Airgun'+str(volumes[ii])])
        TMP = np.fft.fft(np.append(tmp,np.zeros(len(f)-len(tmp))))*np.exp(-ci*2*np.pi*f*(xg[ii]/c))

        x = 2*rtmp*np.cos(angle/180*np.pi) - xg[ii]
        y = yg[ii]
        z = 2*rtmp*np.sin(angle/180*np.pi) - zg[ii]
        rnew = np.sqrt(x**2 + y**2 + z**2)
        p += TMP*np.exp(ci*2*np.pi*f*(rnew/c-tshift))

    if ghost == True:
        for ii in range(len(volumes)):

            tmp = np.copy(AirgunDict['Airgun'+str(volumes[ii])])
            TMP = np.fft.fft(np.append(tmp,np.zeros(len(f)-len(tmp))))*np.exp(-ci*2*np.pi*f*(xg[ii]/c))

            x = 2*rtmp*np.cos(angle/180*np.pi) - xg[ii]
            y = yg[ii]
            z = 2*rtmp*np.sin(angle/180*np.pi) + zg[ii]
            rnew = np.sqrt(x**2 + y**2 + z**2)
            p += -TMP*np.exp(ci*2*np.pi*f*(rnew/c-tshift))

    srcSig = np.fft.ifft(p)

    return srcSig



def CSGtoCMP(data,recs,Shot0=0,Shot1=1400,nShots=97):
    '''
    Rearranges CSG data into CMP data
    INPUTS:
        -data: CSG data for a layer reflection
        -recs: list of receiver channels (same length as data dimension 0)
    OUTPUT:
        -CMPlines: CMP result
    '''

    shots = np.arange(Shot0,Shot1,1)
    CMPshots = np.arange(Shot0,Shot1-nShots,1)

    CMPlines = np.zeros([len(shots)-nShots,len(recs)])

    nRecs = recs[-1]
    cmp_i = 0
    for shot in shots[:-nShots]:
        shotNums = np.linspace(shot-shots[0],shot-shots[0]+nShots,nRecs+1)
        shotNums = shotNums.astype(int)
        CMPlines[cmp_i,:] = data[shotNums[recs],range(len(recs))]
        cmp_i += 1
    
    return CMPlines,CMPshots
    

### PLOT FUNCTIONS ###

def plotSurf(data,x,y,xLims=[0,600],yLims=[0,1400],cLims=[100,200],xlabel='Receiver #',ylabel='Shot #',gauss_sig=0,plotsize=[4,6],title=''):
    

    data[data<cLims[0]] = cLims[0]
    data[data>cLims[1]] = cLims[1]

    data = spnd.gaussian_filter(data,gauss_sig,mode='constant')
    
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    cs = plt.contourf(x,y,data,cmap='turbo',levels=np.linspace(cLims[0],cLims[1],51))
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([xLims[0],xLims[1]])
    plt.ylim([yLims[0],yLims[1]])
    plt.title(title)
    plt.gcf().set_size_inches(plotsize[0],plotsize[1])
    plt.show()
    
    
def plotLine(data,x,shotLims,title='',xlabel='Range (km)',ylabel='Acoustic Energy (dB)',xLims=[0,8],yLims=[100,200]):
    

    dataMean = np.mean(data,axis=0)
    dataStd = np.std(data,axis=0)
    dataMeandB = 20*np.log10(dataMean/(20e-6))
    dataLowdB = 20*np.log10((dataMean-dataStd)/(20e-6))
    dataHighdB = 20*np.log10((dataMean+dataStd)/(20e-6))
    dataScatdB = 20*np.log10(data/(20e-6))
    
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

    for i in range(shotLims[0],shotLims[1]):
        plt.scatter(x,dataScatdB[i,:],c='k',alpha=0.01)
    plt.plot(x,dataMeandB,'r',linewidth=3)
    plt.plot(x,dataLowdB,'r--',linewidth=2)
    plt.plot(x,dataHighdB,'r--',linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xLims)
    plt.ylim(yLims)
    plt.title(title)
    plt.show()
