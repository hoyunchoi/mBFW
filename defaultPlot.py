import matplotlib.pyplot as plt

# Default parameters for plot
#* text
plt.rc('text', usetex = True)
plt.rc('font', size = 50)

#* figure
plt.rc('figure', figsize = (10,10))
plt.rc('axes', linewidth = 3)

#* line
plt.rc('lines', linewidth = 4)
plt.rc('lines', markersize = 15)

#* legend
plt.rc('legend', frameon = False)
plt.rc('legend', fontsize = 40)
plt.rc('legend', handlelength = 1.0)
plt.rc('legend', handletextpad = 0.5)
plt.rc('legend', labelspacing = 0.2)

#* label
plt.rc('axes', labelsize = 50)

#* tick
plt.rc('xtick', direction = 'in')
plt.rc('ytick', direction = 'in')
plt.rc('xtick.major', width = 3)
plt.rc('ytick.major', width = 3)
plt.rc('xtick.major', size = 10)
plt.rc('ytick.major', size = 10)
plt.rc('xtick.major', pad = 10)
plt.rc('ytick.major', pad = 10)
plt.rc('xtick.minor', width = 0.5)
plt.rc('ytick.minor', width = 0.5)

#* tick label
plt.rc('xtick', labelsize = 36)
plt.rc('ytick', labelsize = 36)

#* save
plt.rc('savefig', dpi = 500)
plt.rc('savefig', transparent = True)
plt.rc('savefig', bbox = 'tight')
plt.rc('savefig', format = 'pdf')

if __name__=="__main__":
    print("This is a module draw.py")