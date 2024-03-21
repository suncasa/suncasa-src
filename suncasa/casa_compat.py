"""
casa_compat.py

Provides a uniform interface for importing CASA (Common Astronomy Software Applications) tools across different versions and Python environments. It supports both the monolithic (CASA 4/5/6) and modular (CASA 6+) installations.

This script dynamically imports CASA components, addressing the architectural changes between versions. In monolithic installations, components are available as instantiated objects. In modular installations, components are accessed through `casatools` and `casatasks`.

Function:
- get_casa_tools(alias_list): Returns a dictionary of requested CASA tool instances. It accepts a list of tool aliases and handles dynamic import or built-in object access, ensuring compatibility across CASA versions.

Parameters:
- alias_list (list of str): Aliases for CASA tools to import. Defaults to a common set of tools.

Returns:
- dict: Mapping of tool aliases to their instances or objects.

Usage example:
    casa_tools = get_casa_tools(['tbtool', 'mstool', 'qatool', 'iatool', 'rgtool', 'msmdtool', 'smtool', 'metool'])
    for alias, instance in casa_tools.items():
        print(f"{alias}: {instance}")

The function uses `importlib` for CASA 6+ and falls back to direct access in earlier versions or interactive sessions.
"""

# Attempt dynamic imports based on the version and available modules
try:
    # Use importlib for dynamic imports in CASA 6+ environments
    import importlib
except ImportError:
    # Fallback for environments where importlib is not available (should be rare in practice)
    pass

# Define a mapping between tool aliases and their actual package names in casatools
tool_mapping = {
    'tbtool': 'table',
    'mstool': 'ms',
    'qatool': 'quanta',
    'iatool': 'image',
    'rgtool': 'regionmanager',
    'msmdtool': 'msmetadata',
    'smtool': 'simulator',
    'metool': 'measures',
    'cltool': 'componentlist',
}

def check_dependencies():
    missing_packages = []
    try:
        import aipy
    except ImportError:
        missing_packages.append("aipy-eovsa")

    try:
        import eovsapy
    except ImportError:
        missing_packages.append("eovsapy")

    if missing_packages:
        raise ImportError(
            "The following package(s) are required to use this function: {}. "
            "Please install them to proceed.".format(", ".join(missing_packages))
        )


def import_casatools(alias_list=['tbtool', 'mstool', 'qatool', 'iatool', 'rgtool', 'msmdtool', 'smtool', 'metool']):
    """
    Dynamically imports and returns CASA tools specified by their aliases.

    Parameters:
    alias_list (list of str): Aliases of the CASA tools to be imported and returned.

    Returns:
    dict: A dictionary with keys as tool aliases and values as the imported modules or objects.
    """
    tools = {}
    for alias in alias_list:
        if alias in tool_mapping:
            try:
                try:
                    # Dynamically import the module and create an instance of the tool
                    module = importlib.import_module('casatools')
                    tool_name = tool_mapping[alias]
                    tool_instance = getattr(module, tool_name)
                    tools[alias] = tool_instance
                except:
                    tools[alias] = vars()[alias]
                    # Fallback logic for older Python versions or CASA versions; adjust as needed
                    pass
            except ImportError as e:
                print(f"Error importing {tool_name} as {alias}: {e}")
                # Optionally handle the import error (e.g., skip this tool, log the error, etc.)
        else:
            print(f"No mapping found for alias: {alias}")
            # Optionally handle the case where an alias is not recognized

    return tools

def import_casatasks(*task_names):
    """
    Dynamically imports specified CASA tasks from the casatasks module. This is designed to uniformly import CASA tasks
    for both monolithic and modular CASA installations, addressing the issue where direct task import is not supported
    in monolithic CASA versions.

    Intended for modular CASA 6+ installations where tasks are accessed via casatasks. For monolithic CASA
    (versions 4/5/6), this function provides a unified interface, though direct imports are handled via the global namespace.

    Parameters:
    - task_names (str): Names of CASA tasks to import.

    Returns:
    - dict: Mapping of task names to their corresponding functions.

    Raises:
    - ImportError: If a task is not found in casatasks.
    """

    imported_tasks = {}
    try:
        from casatasks import __dict__ as casatasks_dict
        for task_name in task_names:
            if task_name in casatasks_dict:
                imported_tasks[task_name] = casatasks_dict[task_name]
        # Printing imported tasks
        if imported_tasks:
            print(f"Imported CASA tasks from casatasks: {', '.join(imported_tasks.keys())}.")
        else:
            print("No tasks were imported from casatasks.")
    except ImportError as e:
        print(f"Error importing from casatasks: {e}")
    return imported_tasks