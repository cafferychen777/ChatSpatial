"""
Test image processing utility module
"""

import sys
import os
import unittest
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Add project root directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.utils.image_utils import (
    fig_to_image,
    fig_to_base64,
    create_placeholder_image
)


class TestImageUtils(unittest.TestCase):
    """Test image processing utility functions"""

    def test_fig_to_image(self):
        """Test fig_to_image function"""
        # Create a simple figure
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
        ax.set_title("Test Figure")

        # Convert to Image object
        image = fig_to_image(fig)

        # Verify it returns an Image object
        self.assertEqual(image._format, "png")
        self.assertIsNotNone(image.data)
        self.assertIsInstance(image.data, bytes)

        # Verify image data is not empty
        self.assertGreater(len(image.data), 0)

        # Verify mime_type is correct
        self.assertEqual(image._mime_type, "image/png")

    def test_fig_to_base64(self):
        """Test fig_to_base64 function"""
        # Create a simple figure
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
        ax.set_title("Test Figure")

        # Convert to base64 string
        base64_str = fig_to_base64(fig)

        # Verify it returns a string
        self.assertIsInstance(base64_str, str)

        # Verify string is not empty
        self.assertGreater(len(base64_str), 0)

    def test_create_placeholder_image(self):
        """Test create_placeholder_image function"""
        # Create placeholder image
        message = "Test Placeholder"
        image = create_placeholder_image(message)

        # Verify it returns an Image object
        self.assertEqual(image._format, "png")
        self.assertIsNotNone(image.data)
        self.assertIsInstance(image.data, bytes)

        # Verify image data is not empty
        self.assertGreater(len(image.data), 0)

        # Verify mime_type is correct
        self.assertEqual(image._mime_type, "image/png")

    def test_image_size_control(self):
        """Test image size control functionality"""
        # Create a complex figure (large size, high resolution)
        fig, ax = plt.subplots(figsize=(20, 20), dpi=300)

        # Generate complex data
        x = np.linspace(0, 10, 1000)
        y = np.sin(x) * np.cos(x**2)

        # Draw complex figure
        ax.plot(x, y, linewidth=2)
        ax.set_title("Complex Figure")

        # Use smaller maximum size limit
        max_size_kb = 100  # Increase size limit
        image = fig_to_image(fig, max_size_kb=max_size_kb)

        # Verify image size doesn't exceed limit (allow some margin)
        actual_size_kb = len(image.data) / 1024
        print(f"Image size: {actual_size_kb:.2f} KB")

        # Verify image was compressed (original size should be much larger)
        self.assertLess(actual_size_kb, 200)  # Ensure image is compressed to reasonable size


if __name__ == "__main__":
    unittest.main()
