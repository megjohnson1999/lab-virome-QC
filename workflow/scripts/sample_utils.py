"""
Utility functions for sample detection and management
"""

from pathlib import Path
import sys
import os


def get_input_directories(config):
    """
    Get list of input directories to scan based on configuration.

    Supports three modes:
    1. Single directory (backwards compatible)
    2. Multiple directories
    3. Recursive scanning from a root directory

    Parameters
    ----------
    config : dict
        Pipeline configuration dictionary

    Returns
    -------
    list of Path
        List of directories to scan for samples
    """
    auto_config = config.get("sample_auto_detection", {})

    # Get configuration options
    input_dir = auto_config.get("input_dir")
    input_dirs = auto_config.get("input_dirs")
    recursive = auto_config.get("recursive", False)
    max_depth = auto_config.get("max_depth", None)

    # Validate mutually exclusive options
    if input_dir and input_dirs:
        sys.stderr.write("ERROR: Cannot specify both 'input_dir' and 'input_dirs'. Choose one.\n")
        sys.exit(1)

    if not input_dir and not input_dirs:
        sys.stderr.write("ERROR: Must specify either 'input_dir' or 'input_dirs'\n")
        sys.exit(1)

    if recursive and input_dirs:
        sys.stderr.write("ERROR: Recursive scanning not supported with 'input_dirs'. Use 'input_dir' with 'recursive: true'\n")
        sys.exit(1)

    directories = []

    # Mode 1: Single directory
    if input_dir:
        base_dir = Path(input_dir)

        if not base_dir.exists():
            sys.stderr.write(f"ERROR: Input directory does not exist: {base_dir}\n")
            sys.exit(1)

        if not base_dir.is_dir():
            sys.stderr.write(f"ERROR: Input path is not a directory: {base_dir}\n")
            sys.exit(1)

        if recursive:
            # Recursive scanning
            directories.extend(_get_directories_recursive(base_dir, max_depth))
        else:
            # Single directory (backwards compatible)
            directories.append(base_dir)

    # Mode 2: Multiple directories
    elif input_dirs:
        for dir_path in input_dirs:
            dir_obj = Path(dir_path)

            if not dir_obj.exists():
                sys.stderr.write(f"ERROR: Input directory does not exist: {dir_obj}\n")
                sys.exit(1)

            if not dir_obj.is_dir():
                sys.stderr.write(f"ERROR: Input path is not a directory: {dir_obj}\n")
                sys.exit(1)

            directories.append(dir_obj)

    if not directories:
        sys.stderr.write("ERROR: No valid input directories found\n")
        sys.exit(1)

    print(f"Scanning {len(directories)} director{'y' if len(directories) == 1 else 'ies'} for samples", file=sys.stderr)
    for dir_path in directories:
        print(f"  - {dir_path}", file=sys.stderr)

    return directories


def _get_directories_recursive(root_dir, max_depth=None):
    """
    Recursively find all directories under root_dir.

    Parameters
    ----------
    root_dir : Path
        Root directory to search
    max_depth : int, optional
        Maximum depth to search (None = unlimited)

    Returns
    -------
    list of Path
        List of directories found
    """
    directories = [root_dir]  # Include root directory

    try:
        for item in root_dir.rglob("*"):
            if item.is_dir():
                # Calculate depth relative to root_dir
                relative_path = item.relative_to(root_dir)
                depth = len(relative_path.parts)

                # Skip if exceeds max depth
                if max_depth is not None and depth > max_depth:
                    continue

                directories.append(item)

    except PermissionError as e:
        sys.stderr.write(f"WARNING: Permission denied accessing {e.filename}, skipping...\n")
    except Exception as e:
        sys.stderr.write(f"WARNING: Error accessing directory {root_dir}: {e}\n")

    return sorted(directories)


def scan_directory(input_dir, r1_pattern, r2_pattern):
    """
    Scan a single directory for paired-end samples.

    Parameters
    ----------
    input_dir : Path
        Directory to scan for samples
    r1_pattern : str
        Glob pattern for R1 files (e.g., "*_R1.fastq.gz")
    r2_pattern : str
        Glob pattern for R2 files (e.g., "*_R2.fastq.gz")

    Returns
    -------
    dict
        Dictionary of samples found in this directory
        Format: {sample_name: {"r1": path, "r2": path, "source_dir": path}}
    """
    # Extract the suffix from patterns (everything after the *)
    if "*" not in r1_pattern or "*" not in r2_pattern:
        sys.stderr.write("ERROR: Sample patterns must contain '*' wildcard\n")
        sys.exit(1)

    r1_suffix = r1_pattern.split("*", 1)[1]
    r2_suffix = r2_pattern.split("*", 1)[1]

    # Find all R1 files matching the pattern
    r1_files = sorted(input_dir.glob(r1_pattern))

    if not r1_files:
        return {}  # No warning here, let caller handle empty directories

    samples = {}

    for r1_path in r1_files:
        # Extract sample name by removing the suffix
        if not r1_path.name.endswith(r1_suffix):
            continue

        sample_name = r1_path.name[:-len(r1_suffix)]

        # Look for corresponding R2 file
        r2_path = input_dir / f"{sample_name}{r2_suffix}"

        if r2_path.exists():
            samples[sample_name] = {
                "r1": str(r1_path),
                "r2": str(r2_path),
                "source_dir": str(input_dir)  # Track which directory this came from
            }
            print(f"  Found sample: {sample_name} in {input_dir}", file=sys.stderr)
        else:
            sys.stderr.write(f"WARNING: No R2 file found for sample '{sample_name}' in {input_dir}, skipping...\n")

    return samples


def merge_samples_from_directories(all_directory_samples, conflict_resolution="error"):
    """
    Merge samples from multiple directories, handling naming conflicts.

    Parameters
    ----------
    all_directory_samples : list of dict
        List of sample dictionaries from each directory
    conflict_resolution : str
        Strategy for handling duplicate sample names:
        - "error": Stop with error if duplicates found
        - "prefix_dir": Add directory name prefix
        - "newest": Use sample from most recently modified directory

    Returns
    -------
    dict
        Merged dictionary of samples with conflicts resolved
    """
    merged_samples = {}
    conflicts = {}  # Track conflicts for reporting

    for samples_dict in all_directory_samples:
        for sample_name, sample_info in samples_dict.items():

            if sample_name in merged_samples:
                # Conflict detected!
                if sample_name not in conflicts:
                    conflicts[sample_name] = [merged_samples[sample_name]]
                conflicts[sample_name].append(sample_info)

                if conflict_resolution == "error":
                    # Collect all conflicting directories for error message
                    conflict_dirs = [info["source_dir"] for info in conflicts[sample_name]]
                    sys.stderr.write(f"ERROR: Duplicate sample name '{sample_name}' found in multiple directories:\n")
                    for dir_path in conflict_dirs:
                        sys.stderr.write(f"  - {dir_path}\n")
                    sys.stderr.write("Use conflict_resolution: 'prefix_dir' or 'newest' to handle this automatically\n")
                    sys.exit(1)

                elif conflict_resolution == "newest":
                    # Use the sample from the most recently modified directory
                    existing_sample = merged_samples[sample_name]
                    existing_mtime = os.path.getmtime(existing_sample["source_dir"])
                    new_mtime = os.path.getmtime(sample_info["source_dir"])

                    if new_mtime > existing_mtime:
                        merged_samples[sample_name] = sample_info
                        print(f"  Conflict resolved: Using newer '{sample_name}' from {sample_info['source_dir']}", file=sys.stderr)
                    else:
                        print(f"  Conflict resolved: Keeping existing '{sample_name}' from {existing_sample['source_dir']}", file=sys.stderr)

                elif conflict_resolution == "prefix_dir":
                    # Need to rename both the existing and new sample
                    if sample_name not in [s for s in merged_samples.keys() if "_" in s]:
                        # First time we see this conflict, rename the existing sample too
                        existing_sample = merged_samples[sample_name]
                        existing_dir_name = Path(existing_sample["source_dir"]).name
                        new_existing_name = f"{existing_dir_name}_{sample_name}"

                        # Remove the original and add the prefixed version
                        del merged_samples[sample_name]
                        merged_samples[new_existing_name] = existing_sample
                        print(f"  Conflict resolved: Renamed existing '{sample_name}' to '{new_existing_name}'", file=sys.stderr)

                    # Add the new sample with directory prefix
                    new_dir_name = Path(sample_info["source_dir"]).name
                    new_sample_name = f"{new_dir_name}_{sample_name}"
                    merged_samples[new_sample_name] = sample_info
                    print(f"  Conflict resolved: Renamed new '{sample_name}' to '{new_sample_name}'", file=sys.stderr)

            else:
                # No conflict, add directly
                merged_samples[sample_name] = sample_info

    # Report conflicts summary if any were resolved
    if conflicts and conflict_resolution != "error":
        print(f"\nResolved {len(conflicts)} sample name conflict(s) using '{conflict_resolution}' strategy", file=sys.stderr)

    return merged_samples


def auto_detect_samples(config):
    """
    Auto-detect paired-end samples from directories based on file patterns.

    Now supports multiple directories and recursive scanning.

    Parameters
    ----------
    config : dict
        Pipeline configuration dictionary

    Returns
    -------
    dict
        Dictionary of samples with r1 and r2 paths
        Format: {sample_name: {"r1": path, "r2": path}}
    """
    # If auto-detection is not enabled, return manual samples
    auto_config = config.get("sample_auto_detection", {})
    if not auto_config.get("enabled", False):
        return config.get("samples", {})

    # Get auto-detection settings
    r1_pattern = auto_config.get("r1_pattern", "*_R1.fastq.gz")
    r2_pattern = auto_config.get("r2_pattern", "*_R2.fastq.gz")
    conflict_resolution = auto_config.get("conflict_resolution", "error")

    # Get directories to scan using new helper function
    try:
        directories = get_input_directories(config)
    except SystemExit:
        # Re-raise configuration errors from get_input_directories
        raise

    # Scan each directory for samples
    all_directory_samples = []
    total_directories_with_samples = 0

    for input_dir in directories:
        directory_samples = scan_directory(input_dir, r1_pattern, r2_pattern)

        if directory_samples:
            all_directory_samples.append(directory_samples)
            total_directories_with_samples += 1
            print(f"Found {len(directory_samples)} sample(s) in {input_dir}", file=sys.stderr)
        else:
            print(f"No samples found in {input_dir}", file=sys.stderr)

    # Check if any samples were found
    if not all_directory_samples:
        sys.stderr.write(f"ERROR: No valid paired-end samples found in any of the {len(directories)} director{'y' if len(directories) == 1 else 'ies'}\n")
        sys.stderr.write(f"Search pattern: R1='{r1_pattern}', R2='{r2_pattern}'\n")
        sys.exit(1)

    # Merge samples from all directories, handling conflicts
    samples = merge_samples_from_directories(all_directory_samples, conflict_resolution)

    # Remove source_dir field from final samples (maintain backwards compatibility)
    for sample_name, sample_info in samples.items():
        if "source_dir" in sample_info:
            del sample_info["source_dir"]

    print(f"\nAuto-detected {len(samples)} sample(s) total from {total_directories_with_samples} director{'y' if total_directories_with_samples == 1 else 'ies'}", file=sys.stderr)

    return samples


def get_samples(config):
    """
    Get samples either from manual config or auto-detection.

    This is the main entry point for sample management.

    Parameters
    ----------
    config : dict
        Pipeline configuration dictionary

    Returns
    -------
    dict
        Dictionary of samples with r1 and r2 paths
    """
    return auto_detect_samples(config)
