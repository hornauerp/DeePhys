classdef MEArecording < handle
   properties
       Metadata
       Parameters
       Units
       NetworkFeatures
       UnitFeatures
       Connectivity
       GraphFeatures
   end
   
   methods (Hidden)
       function mearec = MEArecording(metadata, parameters, analyses)
           arguments
              metadata (1,1) struct %Metadata to be parsed // InputPath is required
              parameters (1,1) struct = struct()
              analyses (1,1) struct = struct()
           end
           if nargin > 0
               mearec.Metadata = parseMetadata(metadata);
               
           end
       end
       
       function parseMetadata(obj, metadata)
           if isempty(fieldnames(metadata)) %Check if any metadata is provided
               warning('No metadata provided, cannot continue');
           elseif ~isfield(metadata.InputPath)
               warning('Metadata does not provide InputPath, cannot continue');
               
           else
               md_fields = fieldnames(metadata);
               for f = md_fields
                  obj.(f) = metadata.(f); 
               end
           end
       end
   end
end