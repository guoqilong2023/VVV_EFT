class GetParameter_Range():
    def __init__(self,):
        pass

    def Parameter_Range(self,process="WWW_dim6"):
        if "Dim6" in process:
            return self.dim6_Parameter_Range(process)
        if "Dim8" in process:
            return self.dim8_Parameter_Range(process)

    def dim6_Parameter_Range(self,process):
        Parameter_Range = {} ; Weight_Index = {}
        Parameter_Range['SM'] = [0] ; Weight_Index['SM'] = (0,1)
        Parameter_Range['cW'] = [-10,-5,-1,-0.7,-0.5,-0.3,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.3,0.5,0.7,1,5,10] ; Weight_Index['cW'] = (1,19)
        Parameter_Range['cHbox'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHbox'] = (19,37)
        Parameter_Range['cHDD'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHDD'] = (37,55)
        Parameter_Range['cHW'] = [-10,-5,-1,-0.7,-0.5,-0.3,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.3,0.5,0.7,1,5,10] ; Weight_Index['cHW'] = (55,73)
        # Parameter_Range['cHB'] = [-10,-5,-1,-0.7,-0.5,-0.3,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.3,0.5,0.7,1,5,10] ; Weight_Index['cHB'] = (73,91)
        Parameter_Range['cHB'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHB'] = (73,91)
        Parameter_Range['cHWB'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHWB'] = (91,109)
        Parameter_Range['cHl3'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHl3'] = (109,127)
        Parameter_Range['cHq1'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHq1'] = (127,145)
        Parameter_Range['cHq3'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHq3'] = (145,163)
        Parameter_Range['cHu'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHu'] = (163,181)
        Parameter_Range['cHd'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cHd'] = (181,199)
        Parameter_Range['cll1'] = [-100,-50,-10,-7,-5,-3,-1,-0.5,-0.1,0.1,0.5,1,3,5,7,10,50,100] ; Weight_Index['cll1'] = (199,217)
        return Parameter_Range,Weight_Index

    def dim8_Parameter_Range(self,process):
        if "WWW" in process:
            return self.dim8_WWW_Parameter_Range()
        if "WWZ" in process:
            return self.dim8_WWZ_Parameter_Range()
        if "WZZ" in process:
            return self.dim8_WZZ_Parameter_Range()
        if "ZZZ" in process:
            return self.dim8_ZZZ_Parameter_Range()

    def dim8_WWW_Parameter_Range(self,):
        Parameter_Range = {'FT1': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT0': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT2': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT4': [-3e-06, -1.5e-06, -1e-06, -8e-07, -4e-07, -2e-07, 2e-07, 4e-07, 8e-07, 1e-06, 1.5e-06, 3e-06], 'FM6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM0': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM1': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FS0': [-3e-09, -1.5e-09, -1e-09, -8e-10, -4e-10, -2e-10, 2e-10, 4e-10, 8e-10, 1e-09, 1.5e-09, 3e-09], 'FS1': [-3e-09, -1.5e-09, -1e-09, -8e-10, -4e-10, -2e-10, 2e-10, 4e-10, 8e-10, 1e-09, 1.5e-09, 3e-09], 'FS2': [-3e-09, -1.5e-09, -1e-09, -8e-10, -4e-10, -2e-10, 2e-10, 4e-10, 8e-10, 1e-09, 1.5e-09, 3e-09]}
        Weight_Index = {'FT1': (97, 109), 'FT0': (85, 97), 'FT3': (121, 133), 'FT2': (109, 121), 'FT4': (133, 145), 'FM6': (61, 73), 'FM7': (73, 85), 'FM0': (37, 49), 'FM1': (49, 61), 'FS0': (1, 13), 'FS1': (13, 25), 'FS2': (25, 37)}
        Parameter_Range['SM'] = [0] ; Weight_Index['SM'] = (0,1)
        return Parameter_Range,Weight_Index

    def dim8_WWZ_Parameter_Range(self,):
        Parameter_Range = {'FT1': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT0': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT2': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM2': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM0': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM1': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FS0': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS1': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS2': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08]}
        Weight_Index = {'FT1': (145, 157), 'FT0': (133, 145), 'FT3': (121, 133), 'FT2': (157, 169), 'FT5': (193, 205), 'FT4': (133, 145), 'FT7': (217, 229), 'FT6': (205, 217), 'FM6': (109, 121), 'FM7': (121, 133), 'FM4': (85, 97), 'FM5': (97, 109), 'FM2': (61, 73), 'FM3': (73, 85), 'FM0': (37, 49), 'FM1': (49, 61), 'FS0': (1, 13), 'FS1': (13, 25), 'FS2': (25, 37)}
        Parameter_Range['SM'] = [0] ; Weight_Index['SM'] = (0,1)
        return Parameter_Range,Weight_Index

    def dim8_WZZ_Parameter_Range(self,):
        Parameter_Range = {'FT1': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT0': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT2': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM2': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM0': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM1': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FS0': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS1': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS2': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08]}
        Weight_Index = {'FT1': (145, 157), 'FT0': (133, 145), 'FT3': (121, 133), 'FT2': (157, 169), 'FT5': (193, 205), 'FT4': (133, 145), 'FT7': (217, 229), 'FT6': (205, 217), 'FM6': (109, 121), 'FM7': (121, 133), 'FM4': (85, 97), 'FM5': (97, 109), 'FM2': (61, 73), 'FM3': (73, 85), 'FM0': (37, 49), 'FM1': (49, 61), 'FS0': (1, 13), 'FS1': (13, 25), 'FS2': (25, 37)}
        Parameter_Range['SM'] = [0] ; Weight_Index['SM'] = (0,1)
        return Parameter_Range,Weight_Index

    def dim8_ZZZ_Parameter_Range(self,):
        Parameter_Range = {'FT1': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT0': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT2': [-3e-12, -1.5e-12, -1e-12, -8e-13, -4e-13, -2e-13, 2e-13, 4e-13, 8e-13, 1e-12, 1.5e-12, 3e-12], 'FT5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM6': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM7': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM4': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM5': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM2': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM3': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM0': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FM1': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FS0': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS1': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FS2': [-3e-08, -1.5e-08, -1e-08, -8e-09, -4e-09, -2e-09, 2e-09, 4e-09, 8e-09, 1e-08, 1.5e-08, 3e-08], 'FT8': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11], 'FT9': [-3e-11, -1.5e-11, -1e-11, -8e-12, -4e-12, -2e-12, 2e-12, 4e-12, 8e-12, 1e-11, 1.5e-11, 3e-11]}
        Weight_Index = {'FT1': (145, 157), 'FT0': (133, 145), 'FT3': (121, 133), 'FT2': (157, 169), 'FT5': (193, 205), 'FT4': (133, 145), 'FT7': (217, 229), 'FT6': (205, 217), 'FM6': (109, 121), 'FM7': (121, 133), 'FM4': (85, 97), 'FM5': (97, 109), 'FM2': (61, 73), 'FM3': (73, 85), 'FM0': (37, 49), 'FM1': (49, 61), 'FS0': (1, 13), 'FS1': (13, 25), 'FS2': (25, 37), 'FT8': (229, 241), 'FT9': (241, 253)}
        Parameter_Range['SM'] = [0] ; Weight_Index['SM'] = (0,1)
        return Parameter_Range,Weight_Index

