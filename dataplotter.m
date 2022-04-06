



files = dir('*.txt') ;   
N = length(files) ;

for i = 1:N
    thisfile = files(i).name ;
   

    filename = thisfile
    startRow = 8;

 
    formatSpec = '%*8s%9f%10f%[^\n\r]';

   
    fileID = fopen(filename,'r');


    textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);

   
    fclose(fileID);

 
    CL = dataArray{:, 1};
    CD = dataArray{:, 2};
   
    CD;
    
   
    mikrotero = 10 ;
    megalitero = -10;
    for a = 1 : length(CD);
        if CD(a) < mikrotero;
            mikrotero = CD(a) ;
            imikrotero = a   ;  
            Clideal(i) = CL(a);
            
        end;
        if CL(a) > megalitero ;
            megalitero = CL(a);
            imegalitero = a;
            Clmax(i) = CL(imegalitero);
            
           
            
        end;
        
    end;
    
    

    onoma{i} = files(i).name
    onoma{i} = onoma{i}(1:4)
    
    
    
  
    
    
    


   
    clearvars filename startRow formatSpec fileID dataArray ans;
end
scatter(Clideal,Clmax,'.')
labelpoints(Clideal,Clmax,onoma,'N',0.1)
hold on 



%% Design Point 1  

%airfoil selection 



%inputs 

WprosS = 4.395*4.882 ;        %(N/m^2)
WTOprosS = 4.395*4.882 ;      %(N/m^2)         
Vc = 20 ;                     %(m/s)         
Vs = 12.5 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli1 = Cli  
Clmax_gross1 = Clmax_gross

scatter(Cli1,Clmax_gross1,'*')
labelpoints(Cli1,Clmax_gross1,'DesignPoint1','N',0.1)





%% Design Point 2
%inputs 

WprosS = 2.8*4.882 ;        %(N/m^2)
WTOprosS = 2.8*4.882 ;      %(N/m^2)         
Vc = 20 ;                     %(m/s)         
Vs = 10 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli2 = Cli  
Clmax_gross2 = Clmax_gross

scatter(Cli2,Clmax_gross2,'*')
labelpoints(Cli2,Clmax_gross2,'DesignPoint2','N',0.1)




%% Design Point 3

%inputs 

WprosS = 4.4*4.882 ;        %(N/m^2)
WTOprosS = 4.4*4.882 ;      %(N/m^2)         
Vc = 17.5 ;                     %(m/s)         
Vs = 12.5 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli3 = Cli  
Clmax_gross3 = Clmax_gross

scatter(Cli3,Clmax_gross3,'*')
labelpoints(Cli3,Clmax_gross3,'DesignPoint3','N',0.1)


%% Design Point 4

%inputs 

WprosS = 2.8*4.882 ;        %(N/m^2)
WTOprosS = 2.8*4.882 ;      %(N/m^2)         
Vc = 17.5 ;                     %(m/s)         
Vs = 10 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli4 = Cli  
Clmax_gross4 = Clmax_gross

scatter(Cli4,Clmax_gross4,'*')
labelpoints(Cli4,Clmax_gross4,'DesignPoint4','N',0.1)


%% Design Point 5

%inputs 

WprosS = 4.4*4.882 ;        %(N/m^2)
WTOprosS = 4.4*4.882 ;      %(N/m^2)         
Vc = 15 ;                     %(m/s)         
Vs = 12.5 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli5 = Cli  
Clmax_gross5 = Clmax_gross

scatter(Cli5,Clmax_gross5,'*')
labelpoints(Cli5,Clmax_gross5,'DesignPoint5','N',0.1)


%% Design Point 6

%inputs 

WprosS = 2.8*4.882 ;        %(N/m^2)
WTOprosS = 2.8*4.882 ;      %(N/m^2)         
Vc = 15 ;                     %(m/s)         
Vs = 10 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli6 = Cli  
Clmax_gross6 = Clmax_gross

scatter(Cli6,Clmax_gross6,'*')
labelpoints(Cli6,Clmax_gross6,'DesignPoint6','N',0.1)


%% Design Point 7

%inputs 

WprosS = 1.58*4.882 ;        %(N/m^2)
WTOprosS = 1.58*4.882 ;      %(N/m^2)         
Vc = 15 ;                     %(m/s)         
Vs = 7.5 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli7 = Cli  
Clmax_gross7 = Clmax_gross

scatter(Cli7,Clmax_gross7,'*')
labelpoints(Cli7,Clmax_gross7,'DesignPoint7','N',0.1)


%% Design Point 8

%inputs 

WprosS = 1.58*4.882 ;        %(N/m^2)
WTOprosS = 1.58*4.882 ;      %(N/m^2)         
Vc = 12.5 ;                     %(m/s)         
Vs = 7.5 ;                   %(m/s)      
rgreek = 1.157 ;              %(kg/m^3)   
rgreeksealevel = 1.225 ;      %(kg/m^3)  

%prakseis 

CLc = 2*WprosS/(rgreek*(Vc^2))  ;
CLCw = CLc/0.95  ;
Cli = CLCw/0.9 ;

CLmax = 2*WTOprosS/(rgreeksealevel*(Vs^2)) ;
CLmax_w = CLmax/0.95 ; 
Clmax_gross = CLmax_w/0.9  ;

%outputs

Cli8 = Cli  
Clmax_gross8 = Clmax_gross

scatter(Cli8,Clmax_gross8,'*')
labelpoints(Cli8,Clmax_gross8,'DesignPoint8','N',0.1)

