# From https://stackoverflow.com/questions/11783932/how-to-add-linker-or-compile-flag-in-cmake-file
# This function adds link flags FLAGS to target TARGET
function(target_add_link_flags TARGET FLAGS)
    get_target_property(TEMP ${TARGET} LINK_FLAGS)
    if (TEMP STREQUAL "TEMP-NOTFOUND")
        SET(TEMP "") # set to empty string
    else ()
        SET(TEMP "${TEMP} ") # a space to cleanly separate from existing content
    endif ()
    # append our values
    SET(TEMP "${TEMP}${FLAGS}")
    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS ${TEMP})
endfunction(target_add_link_flags)
