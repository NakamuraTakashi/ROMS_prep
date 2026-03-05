#!/usr/bin/env python
# -*- coding: utf-8 -*-
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import binascii
import math
import time
import gzip

class xrain:
    def __init__(self,input_file,switch=True):
        self.par = self.read_par(input_file)
        self.ppi = self.read_ppi(input_file)
        self.sect = self.read_sect(input_file)
        self.lat = round(self.par[11][0]+(self.par[11][1]/60.0)+(self.par[11][2]/3600.0),4)
        self.lon = round(self.par[12][0]+(self.par[12][1]/60.0)+(self.par[12][2]/3600.0),4)
        self.alt = self.par[13]
        self.range_num = self.par[22]
        self.range_min = self.par[19]
        self.range_max = self.par[20]
        self.range_step = self.par[21]
        self.az_num = self.par[23]
        self.date = self.par[2][0]+"-"+self.par[2][1]+"-"+self.par[2][2]
        self.time = self.par[3][0]+":"+self.par[3][1]+":00"
        self.obs_start = self.par[17][0]+":"+self.par[17][1]+":"+self.par[17][2]
        self.obs_end = self.par[18][0]+":"+self.par[18][1]+":"+self.par[18][2]
        self.elv = self.par[9]
        self.x = self.axis_x()
        self.y = self.axis_y()
        self.z = self.axis_z()
        self.geox = self.axis_geox()
        self.geoy = self.axis_geoy()

        self.rng = self.axis_range()
        if switch==True:
            self.credit()

    def credit(self):
        print("")
        print("  ------------------------------------------------------")
        print("     国土交通省 XRAIN バイナリ解析 Python ライブラリ")
        print("     Version 1.1 (2023年5月5日)")
        print("     作成: 和田有希 (大阪大学大学院工学研究科)")
        print("     https://github.com/YuukiWada/PyXRAIN")
        print("")
        print("     本ソフトウェアを用いて生じた不利益・損害について")
        print("     作成者は一切の責任を負いません。")
        print("  ------------------------------------------------------")
        print("")

    def read_par(self,input_file):
        if ".gz" in input_file:
            f = gzip.open(input_file, "rb")
        else:
            f = open(input_file, "rb")
        data = f.read()
        par = list()
        par.append(binascii.hexlify(data[4:6]).decode())                                    # 0: observation site
        par.append(hex(data[7])[2:])                                                        # 1: data type
        par.append([data[8:12].decode(),data[13:15].decode(),data[16:18].decode()])         # 2: observation date
        par.append([data[19:21].decode(),data[22:24].decode()])                             # 3: observation time
        if data[33]==1:                                                                     # 4: status code
            par.append(0)
        else:
            par.append(1)
        par.append(float(binascii.hexlify(data[40:42]))/10.0)                               # 5: rotation speed
        if binascii.hexlify(data[42:44]).decode()=="0001":                                  # 6: identification of PPI/CAPPI
            par.append("CAPPI")
        else:
            par.append("PPI")
        par.append(int(binascii.hexlify(data[44:46]),16))                                   # 7: elevation step
        par.append(int(binascii.hexlify(data[46:48]),16))                                   # 8: elevation step num
        if int(binascii.hexlify(data[48:50]),16)>32767:                                     # 9: elevation angle
            par.append((int(binascii.hexlify(data[48:50]),16)-2**16)/100.0)
        else:
            par.append(int(binascii.hexlify(data[48:50]),16)/100.0)
        if binascii.hexlify(data[52:56]).decode()=="00000004":                              # 10: site status
            par.append(0)
        else:
            par.append(1)
        par.append([int(binascii.hexlify(data[62:64]),16),int(binascii.hexlify(data[64:66]),16),int(binascii.hexlify(data[66:68]),16)]) # 11: latitude
        par.append([int(binascii.hexlify(data[68:70]),16),int(binascii.hexlify(data[70:72]),16),int(binascii.hexlify(data[72:74]),16)]) # 12: longitude
        par.append(int(binascii.hexlify(data[74:78]),16)/100.0)                             # 13: altitude (m)
        par.append(int(binascii.hexlify(data[88:90]),16)/100.0)                             # 14: horizontal power (kW)
        par.append(int(binascii.hexlify(data[102:104]),16)/100.0)                           # 15: vertical power (kW)
        par.append(int(binascii.hexlify(data[110:112]),16))                                 # 16: transmission frequency (MHz)
        par.append([data[128:130].decode(),data[131:133].decode(),data[134:136].decode()])  # 17: observation start time
        par.append([data[136:138].decode(),data[139:141].decode(),data[142:144].decode()])  # 18: observation start time
        par.append(int(binascii.hexlify(data[144:148]),16)/100.0)                           # 19: range start (m)
        par.append(int(binascii.hexlify(data[148:152]),16)/100.0)                           # 20: range end (m)
        par.append(int(binascii.hexlify(data[152:156]),16)/100.0)                           # 21: range resolution (m)
        par.append(int(binascii.hexlify(data[156:160]),16))                                 # 22: range num
        par.append(int(binascii.hexlify(data[160:162]),16))                                 # 23: azimuth num
        f.close()
        return par

    def read_ppi(self,input_file):
        par = self.par
        if ".gz" in input_file:
            f = gzip.open(input_file, "rb")
        else:
            f = open(input_file, "rb")
        data = f.read()
        num_rng = par[22]
        num_az = par[23]
        value = np.zeros((num_az, num_rng))
        for i in range(num_az):
            for j in range(num_rng):
                index = 512+i*(16+2*num_rng)+(16+2*j)
                value[i,j] = int(binascii.hexlify(data[index:index+2]),16)
        data_type = par[1]
        type_num = ["05", "06", "07", "08", "09", "0E", "11", "12", "15", "19", "21", "25", "31", "35", "9"]
        a = [90.0, 95.0, 100.0, 105.0, 1.0, 80.0, 85.0, 1.0, 1.0, 1.0, 1.0, 1.0, 360.0, 1.0, 1.0]
        b = [0.0, 0.0, 0.0, 0.0, 32768.0, 0.0, 0.0, 32768.0, 32768.0, 1.0, 32768.0, 1.0, 1.0, 32768.0, 32768.0]
        c = [16384, 16384, 16384, 16384, 100.0, 16384, 16384, 100.0, 100.0, 100.0, 100.0, 65533.0, 65534.0, 100.0, 100.0]
        type_index = type_num.index(data_type)
        value = a[type_index]*(value-b[type_index])/c[type_index]
        f.close()
        return value

    def read_sect(self,input_file):
        num_rng = self.par[22]
        num_az = self.par[23]
        sect_info=list()
        if ".gz" in input_file:
            f = gzip.open(input_file, "rb")
        else:
            f = open(input_file, "rb")
        data = f.read()
        for i in range(num_az):
            index = 512+i*(16+2*num_rng)
            sect_info.append([int(binascii.hexlify(data[index:index+2]),16)/100.0,int(binascii.hexlify(data[index+2:index+4]),16)/100.0,\
                              int(binascii.hexlify(data[index+4:index+6]),16)/100.0,int(binascii.hexlify(data[index+6:index+8]),16)/100.0])
        f.close()
        return sect_info

    def axis_x(self):
        n_rng = self.range_num
        n_az = self.az_num
        min_rng = self.range_min
        step_rng = self.range_step
        step_az = 360.0/n_az
        elv = self.elv
        rng = np.tile((min_rng+np.arange(n_rng+1)*step_rng/1000.0)*np.cos(elv*np.pi/180.0), (n_az+1,1))
        az = np.tile(np.arange(n_az+1)*step_az, (n_rng+1,1)).transpose()
        x = rng*np.sin(az*math.pi/180)
        return x

    def axis_y(self):
        n_rng = self.range_num
        n_az = self.az_num
        min_rng = self.range_min
        step_rng = self.range_step
        step_az = 360.0/n_az
        elv = self.elv
        rng = np.tile((min_rng+np.arange(n_rng+1)*step_rng/1000.0)*np.cos(elv*np.pi/180.0), (n_az+1,1))
        #rng = np.tile(min_rng+np.arange(n_rng+1)*step_rng/1000.0, (n_az+1,1))
        az = np.tile(np.arange(n_az+1)*step_az, (n_rng+1,1)).transpose()
        y = rng*np.cos(az*math.pi/180)
        return y

    def axis_geox(self):
        r_earth = 6378.137 # km  
        x = self.x
        lats = self.par[11][0]+self.par[11][1]/60+self.par[11][2]/3600
        lons = self.par[12][0]+self.par[12][1]/60+self.par[12][2]/3600
        km_per_deg = r_earth*np.arccos(np.sin(lats*np.pi/180)*np.sin(lats*np.pi/180)+np.cos(lats*np.pi/180)*np.cos(lats*np.pi/180)*np.cos(np.pi/180))
        geox = lons+x/km_per_deg
        return geox

    def axis_geoy(self):
        r_earth = 6378.137 # km  
        y = self.y
        lats = self.par[11][0]+self.par[11][1]/60+self.par[11][2]/3600
        lons = self.par[12][0]+self.par[12][1]/60+self.par[12][2]/3600
        km_per_deg = r_earth*np.arccos(np.sin(lats*np.pi/180)*np.sin((lats+1)*np.pi/180)+np.cos(lats*np.pi/180)*np.cos((lats+1)*np.pi/180)*np.cos(0))
        geoy = lats+y/km_per_deg
        return geoy
    
    def axis_z(self):
        n_rng = self.range_num
        min_rng = self.range_min
        step_rng = self.range_step
        rng = min_rng+(np.arange(n_rng)+0.5)*step_rng/1000.0
        elv = self.elv
        alt = self.alt
        z = rng*np.tan(elv*math.pi/180.0)+alt/1000.0
        return rng
    
    def axis_range(self):
        n_rng = self.range_num
        min_rng = self.range_min
        step_rng = self.range_step
        rng = min_rng+(np.arange(n_rng)+0.5)*step_rng/1000.0
        return rng
    
    def site(self):
        site_id = ["8105", "8106", "8107", "8108", "8109", "8205", "8206", "8207", "8208", "8209", "820a", "820b",\
                   "8305", "8306", "8405", "8406", "8407", "8408", "8409", "840a", "8505", "8506", "8507", "8508",\
                   "8605", "8606", "8607", "8608", "8609", "860a", "860b", "8705", "8706", "8707", "8708",\
                   "8805", "8806", "8807", "8808", "8900", "8a00", "0000"]

        site_name = ["関東/関東", "関東/新横浜", "関東/氏家", "関東/八丈島", "関東/船橋", "九州/九千地", "九州/菅岳", "九州/古月山", "九州/風師山", "九州/桜島", "九州/山鹿", "九州/宇城",\
                     "北海道/北広島", "北海道/石狩", "東北/一関", "東北/一迫", "東北/涌谷", "東北/岩沼", "東北/伊達", "東北/田村", "北陸/水橋", "北陸/能美", "東北/京ヶ瀬", "東北/中の口",\
                     "中部/尾西", "中部/安城", "中部/鈴鹿", "中部/静岡北", "中部/香貫山", "中部/富士宮", "中部/浜松", "近畿/六甲", "近畿/葛城", "近畿/鷲峰山", "近畿/田口",\
                     "中国/熊山", "中国/常山", "中国/野貝原", "中国/牛尾山", "四国", "沖縄", "未定義"]
        index = site_id.index(self.par[0])
        if index==None:
            return "未定義"
        else:
            return site_name[index]

    def obs_type(self):
        obs_id = ["05", "06", "07", "08", "09", "0E", "11", "12", "15", "19", "21", "25", "31", "35", "9"]
        type_name = ["受信電力強度 [dBm]", "受信電力強度 [dBm]", "受信電力強度 [dBm]", "受信電力強度 [dBm]", "受信電力強度 [dBm]", "受信電力強度 [dBm]", "受信電力強度 [dBm]",
                     "レーダー反射因子 [dBZ]", "風速 [m/s]", "分散 [m/s]", "反射因子差 [dB]", "偏波間相関係数", "偏波間位相差 [deg]", "伝播位相差変化率 [deg/km]",
                     "MTI処理後受信電力強度 [dBm]"]
        index = obs_id.index(self.par[1])
        return type_name[index]

    def parameter(self):
        par = self.par
        status = ["正常","異常"]
        message = list()
        lat = round(par[11][0]+(par[11][1]/60.0)+(par[11][2]/3600.0),4)
        lon = round(par[12][0]+(par[12][1]/60.0)+(par[12][2]/3600.0),4)
        message.append("   レーダーサイト     : "+str(self.site()))
        message.append("   データ種別         : "+str(self.obs_type()))
        message.append("   データ配信日時     : "+str(par[2][0])+"年"+str(par[2][1])+"月"+str(par[2][2])+"日 "+str(par[3][0])+"時"+str(par[3][1])+"分 (日本標準時)")
        message.append("   XRAINステータス    : "+status[par[4]])
        message.append("   サイトステータス   : "+status[par[10]])
        message.append("   レーダー回転速度   : "+str(par[5])+" rpm")
        message.append("   運用モード         : "+par[6])
        message.append("   CAPPI仰角数        : "+str(par[7]))
        message.append("   CAPPI仰角番号      : "+str(par[8]))
        message.append("   CAPPI仰角          : "+str(par[9])+"度")
        message.append("   レーダーサイト緯度 : 北緯 "+str(lat)+"度")
        message.append("   レーダーサイト経度 : 東経"+str(lon)+"度")
        message.append("   レーダーサイト高度 : "+str(par[13])+" m")
        message.append("   水平偏波送信電力   : "+str(par[14])+" kW")
        message.append("   垂直偏波送信電力   : "+str(par[15])+" kW")
        message.append("   送信周波数         : "+str(par[16]/1000)+" GHz")
        message.append("   観測開始時刻       : "+str(par[17][0])+"時"+str(par[17][1])+"分"+str(par[17][2])+"秒 (日本標準時)")
        message.append("   観測終了時刻       : "+str(par[18][0])+"時"+str(par[18][1])+"分"+str(par[18][2])+"秒 (日本標準時)")
        message.append("   視線方向最小距離   : "+str(par[19])+" m")
        message.append("   視線方向最大距離   : "+str(par[20])+" m")
        message.append("   視線方向ステップ   : "+str(par[21])+" m")
        message.append("   視線方向観測数     : "+str(par[22]))
        message.append("   方位角観測数       : "+str(par[23]))
        for string in message:
            print(string)

    def az_start(self,i):
        return self.sect[i][0]

    def az_end(self,i):
        return self.sect[i][1]

    def el_start(self,i):
        return self.sect[i][2]

    def el_end(self,i):
        return self.sect[i][3]

    def plot_color(self):
        colors = ['#FFFFFF', '#A0D2FF', '#218CFF', '#0041FF', '#FAF500', '#FF9900', '#FF2800', '#B40068']
        values = range(len(colors))
        vmax = np.ceil(np.max(values))
        color_list = []
        for v, c in zip(values, colors):
            color_list.append((v/vmax,c))
        custom_color = LinearSegmentedColormap.from_list('custom_cmap', color_list)
        return custom_color
        
class composite:
    def __init__(self,input_file,switch=True):
        if switch==True:
            self.credit()
        self.par = self.read_par(input_file)
        if self.par[0] != 8001:
            print("Error: input file is broken.")
            exit()
        self.date = self.par[1][0]+"-"+self.par[1][1]+"-"+self.par[1][2]
        self.time = self.par[2][0]+":"+self.par[2][1]
        self.nblock = self.par[3]
        self.edge = self.calc_edge(self.par[4],self.par[5])
        self.nmesh = self.calc_nmesh()
        self.mesh = self.calc_mesh()
        self.center = self.calc_center()
        self.comp = self.extract_comp(input_file,self.edge, self.nmesh)

    def credit(self):
        print("")
        print("  ------------------------------------------------------")
        print("     国土交通省 XRAIN 合成雨量データ Python ライブラリ")
        print("     Version 1.0 (2023年5月16日)")
        print("     作成: 和田有希 (大阪大学大学院工学研究科)")
        print("     https://github.com/YuukiWada/PyXRAIN")
        print("")
        print("     本ソフトウェアを用いて生じた不利益・損害について")
        print("     作成者は一切の責任を負いません。")
        print("  ------------------------------------------------------")
        print("")
        
    def read_par(self,input_file):
        if ".gz" in input_file:
            f = gzip.open(input_file, "rb")
        else:
            f = open(input_file, "rb")
        data = f.read()
        par = list()
        par.append(int(binascii.hexlify(data[2:4])))                                        # 0: date type (fixed to #8001)
        par.append([data[8:12].decode(),data[13:15].decode(),data[16:18].decode()])         # 1: observation date
        par.append([data[19:21].decode(),data[22:24].decode()])                             # 2: observation time
        par.append(int(binascii.hexlify(data[42:44]),16))                                   # 3: block number
        par.append(str(int(binascii.hexlify(data[48:50]))))                                # 4: mesh code (sw)
        par.append(str(int(binascii.hexlify(data[50:52]))))                                 # 5: mesh code (ne)
        f.close()
        return par

    def calc_edge(self,sw,ne):
        west = float("1"+sw[2:4])
        east = float("1"+ne[2:4])+1
        south = float(sw[0:2])/1.5
        north = (float(ne[0:2])+1.0)/1.5
        edge = [[west,east],[south,north]]
        return edge

    def calc_nmesh(self):
        edge = self.edge
        nx = round((edge[0][1]-edge[0][0])/(11.25/3600.0))
        ny = round((edge[1][1]-edge[1][0])/(7.5/3600.0))
        nmesh = [nx,ny]
        return nmesh
        
    def calc_mesh(self):
        edge = self.edge
        nmesh = self.nmesh
        x = np.linspace(edge[0][0],edge[0][1],nmesh[0]+1)
        y = np.linspace(edge[1][0],edge[1][1],nmesh[1]+1)
        X,Y = np.meshgrid(x,y)
        return X,Y

    def calc_center(self):
        dx = 11.25/7200.0
        dy = 7.5/7200.0
        edge = self.edge
        nmesh = self.nmesh
        x = np.linspace(edge[0][0]+dx,edge[0][1]-dx,nmesh[0])
        y = np.linspace(edge[1][0]+dy,edge[1][1]-dy,nmesh[1])
        X,Y = np.meshgrid(x,y)
        return X,Y
    
    def extract_comp(self,input_file,edge,nmesh,outside=True):
        comp = np.full((nmesh[1],nmesh[0]),-1.0,dtype=np.float16)
        #comp = np.full((nmesh[1],nmesh[0]),-1.0)
        index = 64
        num_block = self.nblock
        if ".gz" in input_file:
            f = gzip.open(input_file, "rb")
        else:
            f = open(input_file, "rb")
        data = f.read()

        for i in range(num_block):
            code_south = int(binascii.hexlify(data[index:index+1]),16)
            code_west = int(binascii.hexlify(data[index+1:index+2]),16)
            mesh_sec = str(binascii.hexlify(data[index+2:index+3]))
            num_cell = int(binascii.hexlify(data[index+3:index+4]),16)

            delta_x1 = code_west-int(self.par[4][2:4])
            delta_y1 = code_south-int(self.par[4][0:2])
            delta_x2 = int(mesh_sec[3])
            delta_y2 = int(mesh_sec[2])
            start_pos = [320*delta_x1+40*delta_x2, 320*delta_y1+40*delta_y2]
            for j in range(num_cell):
                index_read = index + 4 + 2*j*1600
                rain = np.frombuffer(data[index_read:index_read+3200], dtype=">u2")
                rain = rain.reshape(40,40)
                nx = start_pos[0]+40*j
                ny = start_pos[1]
                rain = np.flipud((rain&0xFFF).astype(np.float16)/10.0)
                comp[ny:ny+40,nx:nx+40] = rain
                #for n in range(40):
                #    for m in range(40):
                #        index_read = index + 4 + 2*j*1600 + 2*n*40 + 2*m
                #        nx = start_pos[0] + 40*j + m
                #        ny = start_pos[1] + (39-n)
                #        comp[ny,nx] = float(int(binascii.hexlify(data[index_read:index_read+2]),16)&0xFFF)/10.0
            index += 4 + 2*num_cell*1600
        np.place(comp, comp>=409.0, -1.0)
        if outside==False:
            np.place(comp, comp<0.0, 0.0)
        return comp

    def lat_index(self,lat):
        dy = 7.5/3600.0
        edge = self.edge[1][0]
        n = int(np.round((lat-edge-dy/2.0)/dy))
        if (n<0) or (n>self.nmesh[1]):
            print("Error: Selected latitude is out of area.")
            print("Erea range: N"+str(self.edge[0][0])+"-"+str(self.edge[0][1])+"deg")
            print("            E"+str(self.edge[1][0])+"-"+str(self.edge[1][1])+"deg")
            exit()
        return n

    def lon_index(self,lon):
        dx = 11.25/3600.0
        edge = self.edge[0][0]
        n = int(np.round((lon-edge-dx/2.0)/dx))
        if (n<0) or (n>self.nmesh[0]):
            print("Error: Selected longitude is out of area.")
            print("Erea range: N"+str(self.edge[0][0])+"-"+str(self.edge[0][1])+"deg")
            print("            E"+str(self.edge[1][0])+"-"+str(self.edge[1][1])+"deg")
            exit()
        return n
