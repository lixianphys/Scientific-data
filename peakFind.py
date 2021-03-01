#Copyright 2021 Lixian WANG. All Rights Reserved.
import pandas as pd

def scissor(c,x,y):
    '''
    Cut the input data into pieces by a critical value c
    Arguments: 
    c: a critical value below which the input data is ignored
    x: input data in x-axis
    y: input data in y-axis. The same length as x
    Return:
    pieces_x: each piece of data in x-axis above c
    pieces_y: each piece of data in y-axis above c
    pieces_id: the index of each piece of data above c in the original input data (list)
    '''
    if len(x) != len(y):
        raise ValueError('x and y should be of the same length')
    else:
        tb = [yy>c for yy in y]
        pieces_x, pieces_y, pieces_id, new_piece_x, new_piece_y, new_piece_id  = [],[],[],[],[],[]
        
        for i, (xx, yy, boolen) in enumerate(zip(x,y,tb)):
                if boolen and i<len(tb)-1: # put into container if boolen is True and not at the end
                    new_piece_x.append(xx)
                    new_piece_y.append(yy)
                    new_piece_id.append(i)
                elif new_piece_x: # if encounter boolen == False and the container is not empty, pack the container
                    pieces_x.append(new_piece_x)
                    pieces_y.append(new_piece_y)
                    pieces_id.append(new_piece_id)
                    new_piece_x, new_piece_y, new_piece_id  = [],[],[] # initiate an empty container
                
    return pieces_x, pieces_y, pieces_id


def peakIden(x,y,pid):
    '''
    Search for peaks in each piece of data
    Arguments:
    x:input data in x-axis.
    y:input data in y-axis. The same length as x
    pid:the index stored
    Return:
    peak_x: x-value of peaks
    peak_y: y-value of peaks
    peak_id: index of peaks
    '''
    if len(x) == 1: # trival case 01
        return x, y, pid
    elif len(x) == 2: # trival case 02
        return [x[[yy==max(y) for yy in y].index(True)]],[max(y)],[pid[[yy==max(y) for yy in y].index(True)]]
    else: # 3-point moving window for three adjacent points
        peak_x, peak_y, peak_id = [],[],[]
        for i,(xx,yy,ppid) in enumerate(zip(x[:-2],y[:-2],pid[:-2])):
            if yy < y[i+1] and y[i+1] > y[i+2]:
                peak_x.append(x[i+1])
                peak_y.append(y[i+1])
                peak_id.append(pid[i+1])
                
        return peak_x, peak_y, peak_id
    

class peakFilter():
    '''
    Call function peakIden to do the job within a range of (xmin,xmax)
    Arguments:
    c: a critical value below which the input data is ignored
    xmin: lower bound of the range
    xmax: upper bound of the range
    Methods:
    peakPos: call peakIden to find peaks within (xmin,xmax)
    markPeak: mark the peak position (x,y) on existing axis handle (ax)
    '''
    def __init__(self,c,xmin,xmax):
        self.c = c
        self.xmin = xmin
        self.xmax = xmax
    def peakPos(self,x,y):
        allpeaks = []
        pieces_x, pieces_y, pieces_id = scissor(self.c,x,y)
        for px, py, pid in zip(pieces_x, pieces_y, pieces_id):
            peak_x, peak_y, peak_id = peakIden(px,py,pid)
            peaks = [{'xv':xv, 'yv':yv, 'pid':pid} for xv, yv,pid in zip(peak_x, peak_y, peak_id) if xv>self.xmin and xv<self.xmax]
            [allpeaks.append(peak) for peak in peaks]
        return pd.DataFrame(allpeaks)
    def markPeak(self,ax,x,y,s,color='k'):
        allpeaks = self.peakPos(x,y)
        if not allpeaks.empty:
            ax.scatter(allpeaks.xv,allpeaks.yv,s=s,marker='v',color=color)
        return None
        