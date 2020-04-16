clear all; 

ds = spreadsheetDatastore('COVID-19-geographic-disbtribution-worldwide.xlsx');
ds.SelectedVariableNames = {'dateRep','cases','geoId','countriesAndTerritories','deaths'};
ds.SelectedVariableTypes(3) = {'categorical'};

data = readall(ds);

%data  = readtable('COVID-19-geographic-disbtribution-worldwide.xlsx');
countrycodes = unique(data.geoId); 
countrydata =[]; 

X = []; 
country = []; 



for iterdmy=1:length(countrycodes)
    
    rows = (data.geoId ==  countrycodes(iterdmy)) ; 
    countrydata = data(rows,:);
    countrydata = sortrows(countrydata,{'dateRep'},{'ascend'}) ; 
    cntname = countrydata.countriesAndTerritories; 
    timeseries = countrydata.cases ; 
    deaths = countrydata.deaths ; 
    dmy_len = length(timeseries); 
    
    if (dmy_len > 40) 
       
    %    figure(iterdmy) ; 
    %    stem(timeseries) ; 
        dmy = zeros(40,1); 
        dmy = timeseries((dmy_len-40+1):dmy_len); 
        dmy = deaths((dmy_len-40+1):dmy_len);
        country = [country;string(cntname(1))] ; 
        X =[X dmy]; 
     %   title(countrycodes(iterdmy)) ; 
    end
    
    
    
end

weights = pca(X'); 

%% compute two features using PCA 
%% each row of matrix Y represents a country 
%% each row contains two features obtained from PCA 
 
Y = X'*weights(:,1:2);  

%Y = [mean(X)' std(X)'] ; 

%Y = tsne(X');
scatter(Y(:,1),Y(:,2));
text(Y(:,1), Y(:,2), country');
xlabel('z1') 
ylabel('z2'); 
