#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple dependency checker for cell communication analysis.

Checks availability of required dependencies and generates a report.
"""

def check_basic_dependencies():
    """Check basic scientific computing dependencies"""
    deps = ['numpy', 'pandas', 'scipy', 'sklearn', 'scanpy']
    results = {}
    
    for dep in deps:
        try:
            __import__(dep)
            results[dep] = True
        except ImportError:
            results[dep] = False
    
    return results

def check_liana():
    """Check LIANA+ availability"""
    try:
        import liana as li
        version = getattr(li, '__version__', 'unknown')
        return True, "LIANA+ version " + str(version) + " available"
    except ImportError as e:
        return False, "LIANA+ not available: " + str(e)

def check_cellphonedb():
    """Check CellPhoneDB availability"""
    try:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
        return True, "CellPhoneDB available"
    except ImportError as e:
        return False, "CellPhoneDB not available: " + str(e)

def generate_report():
    """Generate comprehensive dependency report"""
    print("=" * 60)
    print("CHATSPATIAL CELL COMMUNICATION DEPENDENCY REPORT")
    print("=" * 60)
    
    # Basic dependencies
    print("\nBASIC DEPENDENCIES:")
    basic_deps = check_basic_dependencies()
    all_basic_ok = all(basic_deps.values())
    
    for dep, available in basic_deps.items():
        status = "✅" if available else "❌"
        print("  " + status + " " + dep)
        if not available:
            print("    Install: pip install " + dep)
    
    print("\nBASIC STATUS: " + ("✅ ALL OK" if all_basic_ok else "❌ MISSING DEPENDENCIES"))
    
    # Cell communication methods
    print("\nCELL COMMUNICATION METHODS:")
    
    # LIANA
    liana_ok, liana_msg = check_liana()
    liana_status = "✅" if liana_ok else "❌"
    print("  " + liana_status + " LIANA+ (Primary method)")
    print("    " + liana_msg)
    if not liana_ok:
        print("    Install: pip install liana-py")
    
    # CellPhoneDB
    cellphonedb_ok, cellphonedb_msg = check_cellphonedb()
    cellphonedb_status = "✅" if cellphonedb_ok else "❌"
    print("  " + cellphonedb_status + " CellPhoneDB (Alternative method)")
    print("    " + cellphonedb_msg)
    if not cellphonedb_ok:
        print("    Install: pip install cellphonedb")
    
    # Summary
    methods_available = sum([liana_ok, cellphonedb_ok])
    print("\nOVERALL SUMMARY:")
    
    if all_basic_ok and methods_available > 0:
        print("✅ Ready for cell communication analysis!")
        print("📊 Available methods: " + str(methods_available) + "/2")
    elif all_basic_ok:
        print("⚠️  Basic dependencies OK, but no communication methods available")
        print("   Install at least one: LIANA+ (recommended) or CellPhoneDB")
    else:
        print("❌ Missing critical dependencies")
        print("   Install missing basic dependencies first")
    
    print("\nRECOMMENDED INSTALLATION ORDER:")
    if not all_basic_ok:
        print("1. Install basic dependencies:")
        for dep, available in basic_deps.items():
            if not available:
                print("   pip install " + dep)
    
    if not liana_ok:
        print("2. Install LIANA+ (recommended primary method):")
        print("   pip install liana-py")
    
    if not cellphonedb_ok:
        print("3. Install CellPhoneDB (alternative method):")
        print("   pip install cellphonedb")
    
    print("\n" + "=" * 60)
    
    return {
        'basic_ok': all_basic_ok,
        'liana_ok': liana_ok,
        'cellphonedb_ok': cellphonedb_ok,
        'ready': all_basic_ok and methods_available > 0
    }

if __name__ == "__main__":
    results = generate_report()
    
    # Exit with error code if not ready
    exit_code = 0 if results['ready'] else 1
    exit(exit_code)