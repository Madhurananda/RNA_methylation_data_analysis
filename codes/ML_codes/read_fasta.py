
import sys,os

def read_ref(reffile, interestlist=None, includeseq=True):
   if includeseq:
      seq_list = {}
   else: seq_list = []

   with open( reffile, 'r') as mr:
       t_id = ''
       t_seq = ''
       while True:
           line = mr.readline();
           if not line: break;
           line = line.strip();
           if len(line)==0: continue;

           if line[0] == '>':
               if not t_id=='' and (interestlist==None or t_id in interestlist):
                   if includeseq: seq_list[ t_id ] = t_seq;
                   else: seq_list.append( t_id )
               t_id = line.split()[0];
               if t_id[0] in ['>', '@']: t_id = t_id[1:]
               t_seq = ''
           else:
               if includeseq: t_seq = t_seq + line;
       if len(t_seq)>0:
            if includeseq: seq_list[ t_id ] = t_seq;
            else: seq_list.append( t_id )
   return seq_list


