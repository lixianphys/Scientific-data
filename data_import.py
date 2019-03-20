class sciData_bswp:
            '''define a class for bfield sweep type data'''
            def __init__(self,name,bfield,rxx,rxy,tmp):
                self.name = name
                self.bfield = bfield
                self.rxx = rxx
                self.rxy = rxy
                self.tmp = tmp

class sciData_gswp:
            '''define a closs for gatevoltage sweep type data'''
            def __init__(self,name,gatevoltage,rxx,rxy,tmp):
                self.name = name
                self.gvolt = gatevoltage
                self.rxx = rxx
                self.rxy = rxy
                self.tmp = tmp

def import_col(directory,fname,cols,skiprows,class_type,rxx_sign,rxy_sign):
            '''import multiple files to consititute a new dictionary type data including all
               cols = [bfield, rxx, rxy, temp] for bswp
               cols = [gvolt, rxx, rxy, temp] for gswp'''
            import numpy
            Dic_data = {}
            if class_type == 'bs':
                for f in fname:
                        print(f)
                        read_data = numpy.loadtxt(f, skiprows=skiprows, usecols=cols)
                        read_data = numpy.array(read_data)
                        Dic_data[f] =  sciData_bswp(f, read_data[:, 0], rxx_sign*read_data[:, 1], rxy_sign*read_data[:, 2], read_data[:,3])
            elif class_type =='gs':
                for f in fname:
                        read_data = numpy.loadtxt(f, skiprows=skiprows, usecols=cols)
                        read_data = numpy.array(read_data)
                        Dic_data[f] = sciData_gswp(f, read_data[:, 0], rxx_sign*read_data[:, 1], rxy_sign*read_data[:, 2],read_data[:,3])
            else:
                print('The class_type is not supported in import_col funtion')
                Dic_data = None
            return Dic_data # return the dictionary type data

