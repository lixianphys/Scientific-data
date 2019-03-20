def gswp_plot(Dic_data, fname, Volt_pick, figsize, fig_title):
                '''bswp_plot is defined to plot (Rxx,Rxy)-bfield at a variety of gatevoltages along with Rxx-gatevoltage at zero-field.
                Dic_data = dictionary type data in which each value element is a sciData_bswp class type data
                Volt_pick: picked voltages
                figsize = (float,float) define the size of the entire figure
                fig_title: the title of the figure
                RETURN:
                newf, ax_gateswp, ax_rxx, ax_rxy'''
                import matplotlib.pyplot as plt
                import numpy as np
                # define physical constants
                e = 1.6021766208E-19
                h = 6.62607015E-34
                newf = plt.figure(figsize=figsize) # Create a new figure
                ax_rxx = newf.add_subplot(221)
                ax_rxy = newf.add_subplot(222)
                ax_temp = newf.add_subplot(223)
                ax_tlog = newf.add_subplot(224)

                # plot Rxx-temp dependence
                for volt in Volt_pick:
                       index = 0
                       Temp = [None]*len(fname)
                       R = [None]*len(fname)
                       for f in fname:
                                index_pick = min(range(len(Dic_data[f].gvolt)), key=lambda k: abs(Dic_data[f].gvolt[k]-volt))
                                Temp[index] = np.mean(Dic_data[f].tmp)
                                R[index] = Dic_data[f].rxx[index_pick]
                                index+=1
                       # plot every single curve in outer loop
                       ax_temp.scatter(Temp,R,label='Vg = %02.2fV'% volt)
                       # plot in in Temp**power semilog
                       ax_tlog.scatter([x**-0.5 for x in Temp],R,label='Vg = %02.2fV'% volt)
                       ax_tlog.set_yscale('log')

                ax_temp.set_xlabel('$Temp (K)$')
                ax_temp.set_ylabel('Rxx ($\Omega$)')
                ax_temp.legend()
                ax_tlog.set_xlabel('$Temp^{-1/2} (K)$')
                ax_tlog.set_ylabel('Rxx ($\Omega$)')
                ax_tlog.legend()

                # plot Rxx
                fname = sorted(fname,key = lambda x: np.mean(Dic_data[x].tmp)) # sort fname by its temp value
                for f in fname:
                       ax_rxx.plot(Dic_data[f].gvolt,Dic_data[f].rxx,label='Temp = %02.1fK'% np.mean(Dic_data[f].tmp))
                ax_rxx.set_xlabel('Gate Voltage (V)')
                ax_rxx.set_ylabel('Rxx ($\Omega$)')
                ax_rxx.legend()
                for volt in Volt_pick:
                        ax_rxx.axvline(x=volt, linestyle=':', color='c')

                # plot Rxy
                for f in fname:
                       ax_rxy.plot(Dic_data[f].gvolt,Dic_data[f].rxy,label='Temp = %02.1fK'% np.mean(Dic_data[f].tmp))
                ax_rxy.set_xlabel('Gate Voltage (V)')
                ax_rxy.set_ylabel('Rxy ($\Omega$)')
                ax_rxy.legend()
                for volt in Volt_pick:
                        ax_rxy.axvline(x = volt, linestyle = ':', color = 'c')
                for i in range(2,5):
                        ax_rxy.axhline(y = h/i/e**2,linestyle = ':',color = 'm')
                plt.title(fig_title)
                plt.show()
                return newf, ax_temp, ax_tlog, ax_rxx, ax_rxy