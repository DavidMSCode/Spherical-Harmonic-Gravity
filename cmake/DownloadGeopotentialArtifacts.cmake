if(NOT DEFINED SHG_GEOPOTENTIAL_ARTIFACT_URL)
    message(FATAL_ERROR "SHG_GEOPOTENTIAL_ARTIFACT_URL is not defined")
endif()

if(NOT DEFINED SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE)
    set(SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE "")
endif()

if(NOT DEFINED SHG_GEOPOTENTIAL_ARTIFACTS_DIR)
    message(FATAL_ERROR "SHG_GEOPOTENTIAL_ARTIFACTS_DIR is not defined")
endif()

set(artifact_root "${SHG_GEOPOTENTIAL_ARTIFACTS_DIR}")
set(artifact_dir "${artifact_root}/bin")
set(artifact_archive "${artifact_root}/geopotential_bins.tar.gz")
set(artifact_stamp "${artifact_root}/.download_complete")

file(MAKE_DIRECTORY "${artifact_root}")

if(EXISTS "${artifact_dir}/models.json")
    message(STATUS "Geopotential artifacts already present at ${artifact_dir}")
    file(WRITE "${artifact_stamp}" "present\n")
    return()
endif()

if(SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE AND EXISTS "${SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE}")
    message(STATUS "Using local geopotential artifact bundle at ${SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE}")
    set(artifact_archive "${SHG_GEOPOTENTIAL_ARTIFACT_ARCHIVE}")
else()
    message(STATUS "Downloading geopotential artifact bundle from ${SHG_GEOPOTENTIAL_ARTIFACT_URL}")
    file(DOWNLOAD
        "${SHG_GEOPOTENTIAL_ARTIFACT_URL}"
        "${artifact_archive}"
        SHOW_PROGRESS
        STATUS download_status
        LOG download_log
        TLS_VERIFY ON
    )

    list(GET download_status 0 download_code)
    list(GET download_status 1 download_message)
    if(NOT download_code EQUAL 0)
        message(STATUS "${download_log}")
        message(FATAL_ERROR "Failed to download geopotential artifacts: ${download_message}")
    endif()
endif()

file(REMOVE_RECURSE "${artifact_dir}")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar xzf "${artifact_archive}"
    WORKING_DIRECTORY "${artifact_root}"
    RESULT_VARIABLE extract_result
    OUTPUT_VARIABLE extract_stdout
    ERROR_VARIABLE extract_stderr
)

if(NOT extract_result EQUAL 0)
    message(STATUS "${extract_stdout}")
    message(STATUS "${extract_stderr}")
    message(FATAL_ERROR "Failed to extract geopotential artifact bundle (exit code ${extract_result})")
endif()

if(NOT EXISTS "${artifact_dir}/models.json")
    message(FATAL_ERROR "Geopotential artifact bundle did not produce ${artifact_dir}/models.json")
endif()

file(WRITE "${artifact_stamp}" "downloaded\n")
message(STATUS "Geopotential artifacts ready at ${artifact_dir}")
