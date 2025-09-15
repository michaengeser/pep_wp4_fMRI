 
runs = 6;
subs = 30;
for ir = 1:runs
   for is = 1:subs

       newString = strrep(Params.PublicParams.subjectSource{ir, is}, 'swrsub', 's12wrsub');
       Params.PublicParams.subjectSource{ir, is} = newString;


   end 
end