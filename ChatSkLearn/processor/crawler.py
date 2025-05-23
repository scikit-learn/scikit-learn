# ===================================================================================
# Project: ChatSkLearn
# File: processor/crawler.py
# Description: This file is used to scrape the scikit-learn documentation page and extract the URLs.
# Author: LALAN KUMAR
# Created: [14-04-2025]
# Updated: [29-04-2025]
# LAST MODIFIED BY: LALAN KUMAR [https://github.com/kumar8074]
# Version: 1.0.0
# ===================================================================================
# CAUTION: Please run this file only when there's a update to the documentation page.

import os
import sys
import asyncio
import aiohttp
import re
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urldefrag, urlparse
import json
import random
from aiohttp import TCPConnector

# Dynamically add the project root directory to sys.path
current_file_path = os.path.abspath(__file__)
project_root = os.path.abspath(os.path.join(current_file_path, "../.."))
if project_root not in sys.path:
    sys.path.append(project_root)

MAC_USER_AGENTS = [
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/15.1 Safari/605.1.15",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 13_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 12_4) AppleWebKit/537.36 (KHTML, like Gecko) Firefox/115.0",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 13_2_1) AppleWebKit/605.1.15 (KHTML, like Gecko)",
]

script_dir = os.path.dirname(os.path.abspath(__file__))
print(f"script_dir:{script_dir}")
os.makedirs(script_dir, exist_ok=True)

class DocCrawlerAgent:
    def __init__(self, base_url, max_depth=5, max_pages=10000, concurrency=25):
        self.base_url = base_url.rstrip('/')
        self.base_netloc = urlparse(self.base_url).netloc
        self.base_path_prefix = urlparse(self.base_url).path.strip('/') + '/'
        self.max_depth = max_depth
        self.max_pages = max_pages
        self.concurrency = concurrency  # Number of concurrent requests
        self.visited = set()
        self.failed = {}
        self.to_visit = asyncio.Queue()
        self.results = []
        self.semaphore = asyncio.Semaphore(concurrency)  # Control concurrency
        self.success_file = os.path.join(script_dir, "successful_urls.txt")

        
        # Unwanted patterns to skip
        self.unwanted_patterns = [
            r'/dev/',        # Development version links
            r'/_[^/]+/',     # Internal underscore paths
            r'\.pdf$',       # PDF files
            r'\.zip$',       # ZIP files
            r'/search\.html', # Search page
            r'/_downloads/',  # Downloads
            r'/_images/',     # Images
            r'/_static/',     # Static files
            r'/glossary\.html', # Glossary can be very large
            r'\?highlight='   # Search highlights
        ]
        self.compiled_patterns = [re.compile(pattern) for pattern in self.unwanted_patterns]
        
        # Create or clear the success file at initialization
        with open(self.success_file, 'w') as f:
            f.write("")

    def is_valid_url(self, url):
        parsed = urlparse(url)
        
        # Check if URL matches any unwanted pattern
        for pattern in self.compiled_patterns:
            if pattern.search(url):
                return False
                
        # Ensure it's on the same domain
        is_same_domain = parsed.netloc == self.base_netloc
        
        # Make sure URL contains the base path prefix (e.g., '/stable/')
        path = parsed.path.strip('/')
        contains_prefix = path.startswith(self.base_path_prefix.strip('/')) or not path
        
        return parsed.scheme in {"http", "https"} and is_same_domain and contains_prefix

    def fix_url(self, url):
        """Ensure URLs maintain the proper path structure"""
        parsed = urlparse(url)
        
        if parsed.netloc == self.base_netloc:
            path = parsed.path.strip('/')
            
            # If the path doesn't start with the base prefix (e.g., 'stable')
            if not path.startswith(self.base_path_prefix.strip('/')):
                # Insert the base path prefix
                base_prefix = self.base_path_prefix.strip('/')
                if not path:
                    new_path = f"/{base_prefix}/"
                else:
                    new_path = f"/{base_prefix}/{path}"
                
                # Reconstruct the URL
                return f"{parsed.scheme}://{parsed.netloc}{new_path}{parsed.params}"
        
        return url

    async def fetch(self, session, url):
        try:
            headers = {
                "User-Agent": random.choice(MAC_USER_AGENTS)
            }
            async with self.semaphore:  # Limit concurrent requests
                async with session.get(url, timeout=30, headers=headers, allow_redirects=True) as response:
                    print(f"🌐 Fetch {url} → Status: {response.status}")
                    if response.status == 200:
                        # Append successful URL to the text file
                        with open(self.success_file, 'a') as f:
                            f.write(f"{url}\n")
                        return await response.text()
                    else:
                        self.failed[url] = f"Status {response.status}"
        except Exception as e:
            print(f"❌ Exception fetching {url}: {e}")
            self.failed[url] = str(e)
        return None

    async def worker(self, session, worker_id):
        while len(self.visited) < self.max_pages:
            try:
                url, depth = await asyncio.wait_for(self.to_visit.get(), timeout=5.0)
            except asyncio.TimeoutError:
                # If queue is empty for 5 seconds, check if all workers are idle
                if self.to_visit.empty():
                    break
                continue

            if url in self.visited or depth > self.max_depth:
                self.to_visit.task_done()
                continue

            print(f"🔍 Worker {worker_id} visiting [Depth {depth}]: {url}")
            self.visited.add(url)
            html = await self.fetch(session, url)
            
            if html is None:
                self.to_visit.task_done()
                continue

            self.results.append(url)

            soup = BeautifulSoup(html, "html.parser")
            for link in soup.find_all("a", href=True):
                href = link["href"]
                joined_url = urljoin(url, href)
                clean_url, _ = urldefrag(joined_url)
                
                # Fix URL to maintain proper structure
                fixed_url = self.fix_url(clean_url)
                
                if self.is_valid_url(fixed_url) and fixed_url not in self.visited:
                    await self.to_visit.put((fixed_url, depth + 1))
            
            self.to_visit.task_done()

    async def crawl(self):
        await self.to_visit.put((self.base_url, 0))

        # Optimize connection pooling and DNS cache
        connector = TCPConnector(limit=100, ttl_dns_cache=300)
        async with aiohttp.ClientSession(connector=connector) as session:
            # Create multiple worker tasks
            workers = [asyncio.create_task(self.worker(session, i)) 
                       for i in range(self.concurrency)]
            
            # Wait for all workers to complete
            await asyncio.gather(*workers)

    def save(self, filename="crawl_results.json"):
        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
        with open(filepath, "w") as f:
            json.dump({
                "successful_urls": list(self.results),
                "failed_urls": self.failed,
                "all_attempted": list(self.visited | self.failed.keys())
            }, f, indent=2)

# Usage
if __name__ == "__main__":
    import time
    start_time = time.time()
    
    agent = DocCrawlerAgent("https://scikit-learn.org/stable/", 
                           max_depth=5, 
                           max_pages=10000,
                           concurrency=25)  # Adjust concurrency based on your network and target site
    
    asyncio.run(agent.crawl())
    agent.save("sklearn_crawl_results.json")
    
    duration = time.time() - start_time
    print(f"\n✅ Finished in {duration:.2f} seconds")
    print(f"Crawled {len(agent.results)} URLs (failed: {len(agent.failed)})")
    print(f"Average speed: {len(agent.visited) / duration:.2f} URLs/second")
    print(f"Successfully crawled URLs saved to {agent.success_file}")