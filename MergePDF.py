from PyPDF2 import PdfFileMerger

path1 = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-exc1mv\bsweep\exc1mV_bsweep_ varyUg (02)\Bsweep_at_[-4, -6, -8, 0.3, 0.5, 0.8, -0.3, -0.5]V_exc1mV.pdf'
path2 = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H_exc50mV\bsweep\gate1p2tom0p1 (02)\Bsweep_at_[0.3, 0.6, 0.84, 1.2]V_exc50mV.pdf'
path3 = r'\\system01\QT\Lixian\data\QC 0492 0h\QC 0492 0H-exc500mv\bsweep\gate1p2to0V (02)\Bsweep_at_[0.3, 0.6, 0.84, 0, 1.2]V_exc500mV.pdf'

pdfs = [path1,path2,path3]

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(open(pdf, 'rb'))

with open(r'\\system01\QT\Lixian\data\QC 0492 0h\Exc_1mV_50mV_500mV_cmp.pdf', 'wb') as fout:
    merger.write(fout)