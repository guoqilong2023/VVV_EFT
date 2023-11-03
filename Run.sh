sh Run_SingleUn.sh jes
sh Run_SingleUn.sh jer
sh Run_SingleUn.sh BtagSF
sh Run_SingleUn.sh lepelSF
sh Run_SingleUn.sh lepmuSF
sh Run_SingleUn.sh prefire
sh Run_SingleUn.sh trigger
sh Run_SingleUn.sh pileup

# 计算signal的rate的总系统误差
python signal_UN_Combined.py

# 合并输入文件，同时画图
