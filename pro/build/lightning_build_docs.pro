pro lightning_build_docs
;+
; Name
; ----
;   LIGHTNING_BUILD_DOCS
;
; Purpose
; -------
;   Pulls the procedure and function headers out of the source code, formats
;   them, and drops them in the API reference folder in the lightning
;   documentation. Note that it does not actually 'build' the documentation 
;   with Sphinx. See the notes below for details on how to build with Sphinx.
;
; Calling Sequence
; ----------------
;   ::
;
;       lightning_build_docs
;
; Output
; ------
;   One .rst file per function or procedure, saved in the documentation source
;   directory specified below.
;
; Notes
; -----
;   This procedure itself is not documented in the API reference; this header is
;   provided only as a convenience to the Lightning maintainers.
;
;   The variable ``tree`` in this procedure should be updated to reflect changes
;   to Lightning. Top level functions should be added to the list under the ``top``
;   key in ``tree``. To actually build the documentation with Sphinx installed do::
;
;       cd <path_to_lightning_directory>/docs
;       make html
;       make latexpdf
;
;   This will generate html documentation (``docs/build/html/index.html``) and pdf
;   documentation (``docs/build/latexpdf/``)
;
; Modification History
; --------------------
;   - 2022/03/22: Created (Erik B. Monson)
;   - 2022/06/28: Made lightning_dir a system variable call (Keith Doore)
;   - 2022/06/28: Converted to batch file (Keith Doore)
;   - 2022/06/28: Updated documentation (Keith Doore)
;   - 2022/08/11: Updated with new directory scheme (Keith Doore)
;   - 2022/08/12: Changed from manually having to specify tree to having it search for all files (Keith Doore)
;   - 2022/08/18: Updated sed calls to remove 1-3 leading spaces (Keith Doore)
;-
 On_error, 2
 compile_opt idl2

; Copy pro files into API_reference directory for easy reading
 ; Clear out the API reference directory before copying files
 destination = !lightning_dir + 'docs/source/API_reference/'
 file_delete, destination, /recursive, /quiet
 file_copy, !lightning_dir+'pro', destination, /recursive
 ; Determine file names
 files = file_search(destination, '*.pro')
 ; Remove build and depreciated files as we do not want them in API
 drop_files = where(strmatch(files, '*build*') or strmatch(files, '*depreciated*'), comp=keep_files, /null)
 drop_dirs = file_dirname(files[drop_files])
 drop_dirs = drop_dirs[uniq(drop_dirs, sort(drop_dirs))]
 file_delete, drop_dirs, /recursive, /quiet
 files = files[keep_files]


; Change path to destination so files can be read
 idl_path = pref_get('IDL_PATH')
 pref_set, 'IDL_PATH', '<IDL_DEFAULT>:+/'+destination, /commit


; Create rst files
 foreach file, files do begin
   rst_file_destination = strmid(file, 0, strlen(file)-3)+'rst'
   ; title = f + "$'\n'"
   ; for i = 0, f.strlen() - 1 do title = title + '='
   ; title = title + "$'\n'"
   title = strupcase(file_basename(file, '.pro'))
   title_under = strjoin(replicate('=', strlen(title)))

   ; Perfectly normal, practical usage of sed:
   ; First sed: replace first line with title
   ; Second sed: replace second line with underline for title
   ; Third sed: add a newline after that underline line (could be a problem if a table is in the first column; solution: indent tables)
   ; Fourth sed: if any line starts with one blank space, eat that space
   command = "sed '1s/.*/" + title + "/' doc.tmp | sed '2s/.*/" + title_under +$
             "/' | sed '/^=/G' |" + " sed 's/^ \{1,3\}//g' > " + rst_file_destination

   doc_library, title, print=("cat > doc.tmp" )
   spawn, command
   file_delete, 'doc.tmp'
 endforeach


; Create index.rst files
;   An index file is required for all directories in API_reference
 dirs = strmid(file_dirname(files), strlen(destination))
 split_dirs = strsplit(dirs, '/')
; Determine maximum number of subdirectories in API_reference
 max_depth = 0
 foreach element, split_dirs do max_depth = max([max_depth, n_elements(element)])
 ; Add one extract to max_depth for base level (i.e., destination)
 dir_array = strarr(n_elements(dirs), max_depth+1)
 foreach dir, dirs, i do begin
   split_dir = strsplit(dir, '/', /extract)
   dir_array[i, 1:n_elements(split_dir)] = split_dir
 endforeach

 for i=0, max_depth do begin   
   uniq_dirs = dir_array[uniq(dir_array[*, i], sort(dir_array[*, i])), i]
   if i ne 0 then uniq_dirs = uniq_dirs[where(uniq_dirs ne '', /null)]

   foreach uniq_dir, uniq_dirs, j do begin
     if i ne max_depth then begin
       subdirs = dir_array[where(uniq_dir eq dir_array[*, i]), i+1] 
       subdirs = subdirs[uniq(subdirs, sort(subdirs))]
       subdirs = subdirs[where(subdirs ne '', /null)]
     endif else subdirs = !null
     if i ne 0 then sec_title = uniq_dir else sec_title = 'API Reference'
     sec_title_bar = strjoin(replicate('=', strlen(sec_title)))

     ; Search to see if directory contains pro files
     current_path = destination+strjoin(reform(dir_array[(where(uniq_dir eq dir_array[*, i]))[0], 0:i]),'/')+'/'
     pro_files = file_search(current_path+'*.pro')
     if n_elements(where(pro_files ne '', /null)) ne 0 then has_pro = 1 else has_pro = 0


     openw, index_file, current_path+'/index.rst', /get_lun

     if i eq 0 then begin
       printf, index_file, '.. _api-label:'
       printf, index_file, ''
     endif

     printf, index_file, strupcase(strmid(sec_title, 0, 1))+strmid(sec_title, 1)
     printf, index_file, sec_title_bar
     printf, index_file, ''
     printf, index_file, '.. toctree::'
     printf, index_file, '    :maxdepth: 1'
     printf, index_file, '    :glob:'
     printf, index_file, '    :caption: Contents'
     printf, index_file, ''
     if has_pro then printf, index_file, '    *'
     foreach subdir, subdirs do printf, index_file, '    '+subdir+'/index'
     free_lun, index_file
     close
   endforeach
 endfor


; Move some extra documentation around and reset path.
 file_copy, !lightning_dir + 'filters/README.md', !lightning_dir + 'docs/source/filters.md', /overwrite
 file_delete, files
 pref_set, 'IDL_PATH', idl_path, /commit

end