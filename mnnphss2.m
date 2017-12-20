%just running the newtonhss method to find the optimal alpha under which  
% the minimal cputime is obtained!
% jisuan zuiyoucanshu


fidhss2=fopen('try2.txt','a');
for k=4:1:4                              %for delta=0.001,0.01,or 0.1.
  delta=0.1*k;
  fprintf(fidhss2,'%4.4f\n',delta);
  for q=100:20:100                            % that is for q=10^j! 
     for l=30:10:50                           % that is n!
          [it1,ti1,ot1,t1]=mn_nphss(l,q,2.0,delta);
%fprintf(fidhss2,'%4d  %4d %8.2f %5d %4d %9.4f\n',l,q,9.0,it1,ot1,t1);
          for i=2.1:0.1:3.1 
              [it2,ti2,ot2,t2]=mn_nphss(l,q,i,delta);
%fprintf(fidhss2,'%4d  %4d %8.2f %5d %4d %9.4f\n',l,q,i,it2,ot2,t2);
              if(t2<t1)     
                t1=t2;          %t1 stands for the minimal cputime!
                it1=it2;
                ti1=ti2;
                ot1=ot2;
                alpha=i;
              end
          end
          tt=0;
          for j=1:1:5
              [it,ti,ot,t]=mn_nphss(l,q,alpha,delta); 
              tt=t+tt;
          end
          tt=tt/5;
          fprintf(fidhss2,'%4d  %4d %8.2f %5d %5d %4d %9.4f\n',l,q,alpha,it1,ti1,ot1,tt);
      end
  end
end
fclose(fidhss2);

