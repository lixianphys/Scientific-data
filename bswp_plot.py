def bswp_plot(Dic_data, fname, Volt_pick, add_path, figsize, skiprows, cols,fig_title, rxx_sign):
            '''bswp_plot is defined to plot (Rxx,Rxy)-bfield at a variety of gatevoltages along with Rxx-gatevoltage at zero-field.
            Dic_data = dictionary type data in which each value element is a sciData_bswp class type data
            Volt_pick: picked voltages
            add_path: the path points to gatesweep@zero-field file
            figsize = (float,float) define the size of the entire figure
            skiprows int
            cols = (int1,int2) int1 for the column of gatevoltage, int2 for the column of Rxx
            fig_title: the title of the figure
            rxx_sign: customizable multification factor for rxx
            RETURN:
            newf, ax_gateswp, ax_rxx, ax_rxy'''
            import matplotlib.pyplot as plt
            import numpy as np
            # define physical constants
            e = 1.6021766208E-19
            h = 6.62607015E-34

            newf = plt.figure(figsize=figsize)  # Create a new figure
            ax_rxx = plt.subplot(222)
            ax_rxy = plt.subplot(224)
            ax_gateswp = plt.subplot(121)
            # plot Rxx-Bfield dependence
            index = 0
            for f in fname:
                        ax_rxx.plot(Dic_data[f].bfield, Dic_data[f].rxx, label='$U_g = %02.1fV$' % Volt_pick[index])
                        index += 1
            ax_rxx.set_xlabel('B Field (T)')
            ax_rxx.set_ylabel('Rxx ($\Omega$)')
            ax_rxx.legend()

            # plot Rxy-Bfield dependence
            index = 0
            for f in fname:
                        ax_rxy.plot(Dic_data[f].bfield, Dic_data[f].rxy, label='$U_g = %02.1fV$' % Volt_pick[index])
                        index += 1
            ax_rxy.set_xlabel('B Field (T)')
            ax_rxy.set_ylabel('Rxy ($\Omega$)')
            ax_rxy.legend()
            #
            for i in range(1, 5):
                        ax_rxy.axhline(y=h / i / e ** 2, linestyle=':', color='m')

            # plot Rxx-GateVolt dependence at zero-field conodition
            gateswp_data = np.loadtxt(add_path, skiprows=skiprows, usecols=cols)
            gateswp_data = np.array(gateswp_data)
            ax_gateswp.plot(gateswp_data[:, 0], rxx_sign*gateswp_data[:, 1], color='g')

            for volt in Volt_pick:
                        ax_gateswp.axvline(x=volt, linestyle=':', color='c')
            ax_gateswp.set_xlabel('$U_g (V)$')
            ax_gateswp.set_ylabel('Rxx ($\Omega$)')

            # Scatter Rxx(at zero-field)
            R = [None] * len(Volt_pick)
            index = 0
            for f in fname:
                        index_pick = min(range(len(Dic_data[f].bfield)), key=lambda k: abs(Dic_data[f].bfield[k]))
                        R[index] = Dic_data[f].rxx[index_pick]
                        index += 1
            # plot every single curve in outer loop
            ax_gateswp.scatter(Volt_pick, R)
            plt.title(fig_title)
            plt.show()
            return newf,ax_gateswp, ax_rxx, ax_rxy