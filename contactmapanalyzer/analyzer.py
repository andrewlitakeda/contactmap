import os
import io
import tempfile
import pandas as pd
from Bio.PDB import PDBParser
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

AA_DICT = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'HSE': 'H', 'HSD': 'H', 'HSP': 'H', 'SEC': 'U',
    'PYL': 'O', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'XLE': 'J',
    'UNK': 'X'
}

class ContactMapAnalyzer:
    def __init__(self, cutoff=4.0):
        self.cutoff = cutoff
        self.parser = PDBParser(QUIET=True)

    def load_structure(self, pdb_input, structure_id="structure"):
        """
        Loads a PDB structure from a file path or file-like object.

        Args:
            pdb_input (str or file-like object): Path to the PDB file or a file-like object.
            structure_id (str): The ID to assign to the structure.

        Returns:
            Bio.PDB.Structure.Structure: The loaded PDB structure object.

        Raises:
            FileNotFoundError: If pdb_input is a path and the file does not exist.
            Exception: For other PDB parsing errors.
        """
        try:
            # If pdb_input is a string path, check if it exists
            if isinstance(pdb_input, str):
                if not os.path.exists(pdb_input):
                    raise FileNotFoundError(f"PDB file not found at: {pdb_input}")
            # PDBParser can handle both file paths (strings) and file-like objects
            structure = self.parser.get_structure(structure_id, pdb_input)
            logger.info(f"Successfully loaded structure '{structure_id}'")
            return structure
        except FileNotFoundError as fnf_error:
             logger.error(f"PDB file loading error: {fnf_error}")
             raise # Re-raise the specific error
        except Exception as e:
            logger.error(f"Failed to load structure '{structure_id}': {e}")
            raise # Re-raise other exceptions

    def get_chain_residues(self, structure, chain_id):
        """
        Extracts residues, their IDs, and names for a specific chain from a structure.

        Args:
            structure (Bio.PDB.Structure.Structure): The PDB structure object.
            chain_id (str): The ID of the chain to process.

        Returns:
            tuple: A tuple containing:
                - list: List of Bio.PDB.Residue.Residue objects.
                - list: List of residue sequence identifiers (integers).
                - list: List of residue names (e.g., 'ALA', 'GLY').

        Raises:
            ValueError: If the specified chain_id is not found in the structure's first model.
        """
        if not structure or not structure[0].has_id(chain_id):
             raise ValueError(f"Chain '{chain_id}' not found in the first model of the structure.")
        residues, residue_ids, residue_names = [], [], []
        for residue in structure[0][chain_id]:
            # Standard residues have hetfield ' '
            # The residue id tuple is (hetfield, sequence_identifier, insertion_code)
            if residue.id[0] == ' ':
                residues.append(residue)
                residue_ids.append(residue.id[1]) # sequence identifier
                residue_names.append(residue.resname)
        logger.info(f"Found {len(residues)} standard residues in chain '{chain_id}'.")
        return residues, residue_ids, residue_names

    def calculate_contact_map(self, chain1_residues, chain2_residues, chain1_ids, chain2_ids, chain1_names, chain2_names):
        """
        Calculates contacts between two sets of residues based on minimum atom distance.

        Args:
            chain1_residues (list): List of Residue objects for chain 1.
            chain2_residues (list): List of Residue objects for chain 2.
            chain1_ids (list): List of residue IDs for chain 1.
            chain2_ids (list): List of residue IDs for chain 2.
            chain1_names (list): List of residue names for chain 1.
            chain2_names (list): List of residue names for chain 2.

        Returns:
            tuple: A tuple containing:
                - dict: A dictionary representing the contact map (e.g., {res1_id: {res2_id: 1}}).
                - list: A list of tuples, each representing a contact:
                        (res1_id, res2_id, res1_amino_acid_code, res2_amino_acid_code).
                        Returns an empty list if no contacts are found or input chains are empty.
        """
        contact_list = []
        if not chain1_residues or not chain2_residues:
            logger.warning("One or both residue lists are empty. No contacts can be calculated.")
            return {}, []

        num_comparisons = len(chain1_residues) * len(chain2_residues)
        logger.info(f"Calculating contacts based on cutoff {self.cutoff} Å for {num_comparisons} residue pairs.")

        for i, res1 in enumerate(chain1_residues):
            res1_id, res1_name = chain1_ids[i], chain1_names[i]
            res1_code = AA_DICT.get(res1_name, 'X')
            for j, res2 in enumerate(chain2_residues):
                res2_id, res2_name = chain2_ids[j], chain2_names[j]
                res2_code = AA_DICT.get(res2_name, 'X')
                min_distance = float('inf')
                try:
                    # Find the minimum distance between any atom pair in the two residues
                    for atom1 in res1:
                        for atom2 in res2:
                            # Bio.PDB atom subtraction calculates distance
                            distance = atom1 - atom2
                            if distance < min_distance:
                                min_distance = distance
                except Exception as dist_err:
                    # Log specific distance calculation errors but continue
                    logger.warning(f"Could not calculate distance between residue {res1_id} (chain1) and {res2_id} (chain2): {dist_err}")
                    continue # Skip this residue pair

                # Check if the minimum distance is within the cutoff
                if min_distance <= self.cutoff:
                    contact_list.append((res1_id, res2_id, res1_code, res2_code))

        # Create a dictionary format for the contact map as well
        contact_dict = {}
        for r1, r2, _, _ in contact_list:
             contact_dict.setdefault(r1, {})[r2] = 1 # Use 1 to indicate contact

        logger.info(f"Found {len(contact_list)} contacts.")
        return contact_dict, contact_list

    def export_to_excel_stream(self, results, output_stream):
        """
        Exports analysis results (contacts per chain pair) to an Excel file stream.

        Args:
            results (dict): A dictionary where keys are interface names (e.g., "A-B")
                            and values are dictionaries containing 'contact_list',
                            'n_residues1', 'n_residues2'.
            output_stream (io.BytesIO or similar): A writable binary stream to write the Excel file to.

        Raises:
            ImportError: If 'xlsxwriter' is not installed.
            Exception: For errors during Excel writing.
        """
        try:
            import xlsxwriter # Ensure xlsxwriter is available
        except ImportError:
             logger.error("The 'xlsxwriter' package is required to export to Excel. Please install it.")
             raise

        logger.info(f"Exporting {len(results)} interface results to Excel stream.")
        try:
            with pd.ExcelWriter(output_stream, engine='xlsxwriter') as writer:
                workbook = writer.book
                summary_data = {"Chain Interface": [], "Residues Chain 1": [], "Residues Chain 2": [], "Contacts Found": []}

                # Define formats
                header_format = workbook.add_format({'bold': True, 'bg_color': '#D7E4BC', 'border': 1, 'align': 'center', 'valign': 'vcenter'})
                cell_format = workbook.add_format({'border': 1, 'align': 'center'}) # Basic cell format

                # Process each interface
                for interface_name, data in results.items():
                    chain1, chain2 = interface_name.split('-') # Assumes "C1-C2" format
                    contact_list = data.get("contact_list", [])
                    n_contacts = len(contact_list)

                    # Add to summary
                    summary_data["Chain Interface"].append(interface_name)
                    summary_data["Residues Chain 1"].append(data.get("n_residues1", "N/A"))
                    summary_data["Residues Chain 2"].append(data.get("n_residues2", "N/A"))
                    summary_data["Contacts Found"].append(n_contacts)

                    if n_contacts == 0:
                        logger.info(f"No contacts found for interface {interface_name}, skipping detailed sheet.")
                        continue # Don't create a sheet if there are no contacts

                    # Prepare DataFrame for the detailed sheet
                    df = pd.DataFrame(
                        contact_list,
                        columns=[f"Residue Chain {chain1}", f"Residue Chain {chain2}", f"AA Chain {chain1}", f"AA Chain {chain2}"]
                    )[[f"Residue Chain {chain1}", f"AA Chain {chain1}", f"Residue Chain {chain2}", f"AA Chain {chain2}"]] # Ensure column order

                    # Write DataFrame to a sheet named after the interface
                    df.to_excel(writer, sheet_name=interface_name, index=False, startrow=1, header=False) # Start writing data from row 1
                    sheet = writer.sheets[interface_name]

                    # Write header row with formatting and set column widths
                    for col_num, value in enumerate(df.columns.values):
                        # Estimate column width needed
                        header_len = len(str(value))
                        # Calculate max length in the column, handle potential non-string types and empty df
                        try:
                             max_len_data = df[value].astype(str).map(len).max() if not df.empty else 0
                        except Exception: # Fallback if astype/map fails
                            max_len_data = 10
                        col_width = max(header_len, max_len_data, 10) + 2 # Add padding
                        sheet.set_column(col_num, col_num, col_width, cell_format) # Apply basic format to column
                        sheet.write(0, col_num, value, header_format) # Write header


                # Create Summary Sheet
                if not summary_data["Chain Interface"]:
                     logger.warning("No interfaces processed, summary sheet will be empty.")

                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name="Summary", index=False, startrow=1, header=False) # Start writing data from row 1
                summary_sheet = writer.sheets["Summary"]

                # Write header row for summary and set column widths
                for col_num, value in enumerate(summary_df.columns.values):
                     header_len = len(str(value))
                     try:
                        max_len_data = summary_df[value].astype(str).map(len).max() if not summary_df.empty else 0
                     except Exception:
                        max_len_data = 10
                     col_width = max(header_len, max_len_data, 15) + 2 # Add padding
                     summary_sheet.set_column(col_num, col_num, col_width, cell_format) # Apply basic format to column
                     summary_sheet.write(0, col_num, value, header_format) # Write header

                logger.info("Excel export process completed.")
        except Exception as e:
            logger.error(f"Failed to export results to Excel: {e}")
            raise # Re-raise exception after logging

def parse_residue_selection(selection_str, all_residue_ids):
    """
    Parses a string to select residue IDs from a list of available IDs.

    The selection string can contain comma-separated numbers or ranges (e.g., "1-10, 25, 50-").
    '-' at the start means from the minimum available ID.
    '-' at the end means up to the maximum available ID.

    Args:
        selection_str (str): The string containing residue selections. If empty or None,
                             returns all available residue IDs.
        all_residue_ids (list or set): A collection of available integer residue IDs.

    Returns:
        set: A set of selected integer residue IDs.

    Raises:
        ValueError: If the selection string format is invalid or refers to IDs not
                    present in all_residue_ids (unless the string is empty).
        TypeError: If all_residue_ids contains non-integer elements.
    """
    if not isinstance(all_residue_ids, (list, set)):
        raise TypeError("all_residue_ids must be a list or set.")

    # Ensure all residue IDs are integers and handle empty case
    try:
        all_residue_ids_set = set(int(r_id) for r_id in all_residue_ids)
    except (ValueError, TypeError) as e:
         raise TypeError(f"all_residue_ids must contain only integer residue numbers. Error: {e}")

    if not selection_str or selection_str.strip() == "":
        logger.debug("No residue selection string provided, returning all available IDs.")
        return all_residue_ids_set # Return all if selection is empty

    if not all_residue_ids_set:
        if selection_str.strip():
             raise ValueError("Residue selection string provided, but no standard residues are available in the input.")
        return set() # Return empty set if no residues and no selection

    selected_ids = set()
    parts = re.split(r'\s*,\s*', selection_str.strip()) # Split by comma, allowing whitespace

    min_res_id = min(all_residue_ids_set)
    max_res_id = max(all_residue_ids_set)
    logger.debug(f"Parsing selection string: '{selection_str}'. Available residue range: {min_res_id}-{max_res_id}")

    for part in parts:
        part = part.strip()
        if not part: continue # Skip empty parts resulting from trailing commas etc.

        # Regex to match ranges like "1-10", "-10", "50-", or single numbers "25"
        range_match = re.match(r'^(-?\d+)?-(-?\d+)?$', part)
        single_match = re.match(r'^(-?\d+)$', part)

        if range_match:
            start_str, end_str = range_match.groups()

            # If only a hyphen is given, it's invalid
            if start_str is None and end_str is None:
                 raise ValueError(f"Invalid range format: '{part}'. Hyphen requires at least one number.")

            # Determine start and end of the range, defaulting to min/max
            try:
                start = int(start_str) if start_str is not None else min_res_id
                end = int(end_str) if end_str is not None else max_res_id
            except ValueError:
                raise ValueError(f"Invalid number in range: '{part}'")

            if start > end:
                raise ValueError(f"Invalid range: start ({start}) is greater than end ({end}) in '{part}'")

            # Add all residues within the specified range that are also in the available set
            found_in_range = {r for r in all_residue_ids_set if start <= r <= end}
            if not found_in_range:
                 logger.warning(f"Selection range '{part}' (interpreted as {start}-{end}) did not match any available residues.")
            selected_ids.update(found_in_range)
            logger.debug(f"Range '{part}' selected {len(found_in_range)} residues: {sorted(list(found_in_range))[:10]}...") # Log first few

        elif single_match:
            try:
                res_id = int(single_match.group(1))
                if res_id in all_residue_ids_set:
                    selected_ids.add(res_id)
                    logger.debug(f"Selected single residue: {res_id}")
                else:
                     logger.warning(f"Selected residue ID {res_id} is not present in the available residues.")
            except ValueError:
                # This should theoretically not happen if regex matched, but good practice
                raise ValueError(f"Invalid number format: '{part}'")
        else:
            raise ValueError(f"Invalid format for selection part: '{part}'. Use numbers or ranges like '1-10'.")

    # Final check: if a non-empty selection string resulted in zero selected IDs
    if not selected_ids and selection_str.strip():
        # This case might occur if the selection contained only IDs not present
        logger.warning(f"The residue selection string '{selection_str}' did not match any of the available residues.")
        # Depending on desired behavior, you might want to raise an error here or return empty set.
        # Current behavior: return empty set, warning already logged.

    logger.info(f"Residue selection parsed. Selected {len(selected_ids)} residues.")
    return selected_ids

def run_contact_analysis(pdb_path, cutoff=4.0, output_filename=None):
    """
    Performs a full contact map analysis for all unique chain pairs in a PDB file
    and saves the results to an Excel file.

    Args:
        pdb_path (str): The full path to the input PDB file.
        cutoff (float, optional): The distance cutoff in Angstroms for defining a contact. Defaults to 4.0.
        output_filename (str, optional): The desired name for the output Excel file.
            If None, it defaults to "<pdb_filename_base>_contacts_<cutoff>A.xlsx"
            and is saved in the same directory as the input PDB file.

    Returns:
        str: The full path to the generated Excel report file, or None if analysis failed.

    Raises:
        FileNotFoundError: If the pdb_path does not exist.
        ValueError: If the PDB file contains fewer than two chains.
        Exception: Catches and logs other potential errors during processing.
    """
    logger.info(f"Starting contact analysis for PDB: {pdb_path} with cutoff: {cutoff} Å")

    analyzer = ContactMapAnalyzer(cutoff=cutoff)

    try:
        # --- 1. Load Structure ---
        if not os.path.exists(pdb_path):
             raise FileNotFoundError(f"Input PDB file not found: {pdb_path}")
        pdb_basename = os.path.basename(pdb_path)
        pdb_name_only, _ = os.path.splitext(pdb_basename)
        structure_id = pdb_name_only # Use PDB filename base as structure ID
        structure = analyzer.load_structure(pdb_path, structure_id=structure_id)

        # --- 2. Identify Chains ---
        available_chains = sorted([c.id for c in structure[0].get_chains() if any(r.id[0] == ' ' for r in c)]) # Only include chains with standard residues
        if len(available_chains) < 2:
            logger.warning(f"PDB file '{pdb_basename}' contains fewer than two chains with standard residues ({available_chains}). Cannot perform inter-chain analysis.")
            return None # Or raise ValueError("Need at least two chains for inter-chain analysis.")
        logger.info(f"Found chains with standard residues: {', '.join(available_chains)}")

        # --- 3. Generate Chain Pairs ---
        chain_pairs = []
        for i, c1 in enumerate(available_chains):
            for c2 in available_chains[i+1:]:
                chain_pairs.append(tuple(sorted((c1, c2)))) # Ensure consistent pair order (e.g., ('A', 'B'))
        logger.info(f"Analyzing {len(chain_pairs)} chain pairs: {', '.join([f'{p[0]}-{p[1]}' for p in chain_pairs])}")

        # --- 4. Calculate Contacts for Each Pair ---
        results = {}
        all_chain_data = {} # Cache residue data

        for chain1_id, chain2_id in chain_pairs:
            interface_name = f"{chain1_id}-{chain2_id}"
            logger.info(f"Processing interface: {interface_name}")

            try:
                # Get residues (use cache)
                if chain1_id not in all_chain_data:
                    all_chain_data[chain1_id] = analyzer.get_chain_residues(structure, chain1_id)
                if chain2_id not in all_chain_data:
                    all_chain_data[chain2_id] = analyzer.get_chain_residues(structure, chain2_id)

                chain1_residues, chain1_ids, chain1_names = all_chain_data[chain1_id]
                chain2_residues, chain2_ids, chain2_names = all_chain_data[chain2_id]

                # Calculate contacts
                contact_dict, contact_list = analyzer.calculate_contact_map(
                    chain1_residues, chain2_residues,
                    chain1_ids, chain2_ids,
                    chain1_names, chain2_names
                )

                results[interface_name] = {
                    "contact_list": contact_list,
                    "n_residues1": len(chain1_ids),
                    "n_residues2": len(chain2_ids),
                }
                logger.info(f"Found {len(contact_list)} contacts for {interface_name}.")

            except ValueError as ve:
                 logger.error(f"Skipping pair {interface_name} due to error getting residues: {ve}")
                 continue # Skip pair if residue fetching fails
            except Exception as e:
                 logger.error(f"Skipping pair {interface_name} due to unexpected error during calculation: {e}")
                 continue # Skip pair on other calculation errors

        # --- 5. Determine Output Path ---
        output_dir = os.path.dirname(pdb_path)
        if output_filename is None:
            # Default filename: <pdb_name>_contacts_<cutoff>A.xlsx
            safe_cutoff_str = str(cutoff).replace('.', 'p') # Make cutoff safe for filename
            default_filename = f"{pdb_name_only}_contacts_{safe_cutoff_str}A.xlsx"
            output_path = os.path.join(output_dir, default_filename)
        else:
            # Use provided filename, ensuring it's in the same directory
            output_path = os.path.join(output_dir, os.path.basename(output_filename))

        # --- 6. Export to Excel ---
        if not results:
             logger.warning("No results were generated (perhaps no valid pairs or contacts found). Skipping Excel export.")
             return None

        excel_stream = io.BytesIO()
        try:
            analyzer.export_to_excel_stream(results, excel_stream)
            excel_stream.seek(0) # Rewind stream

            # Save the stream to file
            with open(output_path, "wb") as f:
                f.write(excel_stream.read())
            logger.info(f"Contact map report successfully saved to: {output_path}")
            return output_path # Return the path of the saved file

        except ImportError as imp_err:
            logger.error(f"Excel export failed: {imp_err}. Make sure 'xlsxwriter' is installed.")
            raise # Re-raise import error
        except Exception as e:
            logger.error(f"Failed to write Excel file to {output_path}: {e}")
            raise # Re-raise other export/write errors

    except FileNotFoundError as e:
         logger.error(f"Analysis failed: {e}")
         raise
    except ValueError as e: # Catch specific errors like chain not found or < 2 chains
         logger.error(f"Analysis prerequisites not met: {e}")
         # Decide whether to raise or return None based on severity
         # For demo, maybe just log and return None is okay if < 2 chains
         # But raise if chain ID logic failed unexpectedly
         if "fewer than two chains" in str(e):
             return None
         else:
             raise
    except Exception as e:
        logger.exception(f"An unexpected error occurred during contact analysis for {pdb_path}: {e}")
        # Depending on desired robustness, might just return None or re-raise
        raise # Re-raise unexpected errors

    return None # Should only be reached if logic error allows skipping return statements above
