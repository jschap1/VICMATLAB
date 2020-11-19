function resultsstruct = make_results_struct_wb(results)

resultsstruct = struct();
resultsstruct.evap = results(:,4);
resultsstruct.runoff = results(:,5);
resultsstruct.baseflow = results(:,6);
resultsstruct.canop_evap = results(:,7);
resultsstruct.transp = results(:,8);
resultsstruct.bare_soil_evap = results(:,9);
resultsstruct.sub_canop = results(:,10);
resultsstruct.sub_snow = results(:,11);
resultsstruct.precip = results(:,13);
resultsstruct.sm0 = results(:,21);
resultsstruct.sm1 = results(:,22);
resultsstruct.sm2 = results(:,23);
resultsstruct.swe = results(:,27);

return